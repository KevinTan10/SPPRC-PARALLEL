#include <boost/config.hpp>
#include <random>
#include <vector>
#include <utility>
#include <chrono>

#ifdef BOOST_MSVC
#pragma warning(disable : 4267)
#endif

#include <boost/graph/adjacency_list.hpp>
#include "r_c_shortest_paths_parallel.hpp"

#include <iostream>
#include <fstream>
#include <immintrin.h>

#define _REAL_DATA_
#define _USE_SIMD_


struct SPPRC_Example_Graph_Vert_Prop
{
    SPPRC_Example_Graph_Vert_Prop(int n = 0)
        : num(n)
    {
    }
    int num;
};

struct SPPRC_Example_Graph_Arc_Prop
{
    SPPRC_Example_Graph_Arc_Prop(int n, std::vector<float> label)
        : num(n), label(label)
    {
    }
    int num;
    std::vector<float> label;
};

typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, SPPRC_Example_Graph_Vert_Prop,
    SPPRC_Example_Graph_Arc_Prop >
    SPPRC_Example_Graph;


// ResourceContainer model
struct spp_spprc_res_cont
{
    spp_spprc_res_cont(std::vector<float> label) : label(label)
    {
    }
    spp_spprc_res_cont& operator=(const spp_spprc_res_cont& other)
    {
        if (this == &other)
            return *this;
        this->~spp_spprc_res_cont();
        new (this) spp_spprc_res_cont(other);
        return *this;
    }
    std::vector<float> label;
};


class ref_spprc
{
public:
    inline void operator()(const SPPRC_Example_Graph& g,
        spp_spprc_res_cont& new_cont, const spp_spprc_res_cont& old_cont,
        boost::graph_traits< SPPRC_Example_Graph >::edge_descriptor ed) const
    {
        const SPPRC_Example_Graph_Arc_Prop& arc_prop = get(boost::edge_bundle, g)[ed];
        const SPPRC_Example_Graph_Vert_Prop& vert_prop
            = get(boost::vertex_bundle, g)[target(ed, g)];
#ifndef _USE_SIMD_
        for (int i = 0; i < new_cont.label.size(); ++i) {
            new_cont.label[i] = old_cont.label[i] + arc_prop.label[i];
        }
#else
        for (int i = 0; i < new_cont.label.size(); i += 8) {
            __m256 v1 = _mm256_loadu_ps(&old_cont.label[i]);
            __m256 v2 = _mm256_loadu_ps(&arc_prop.label[i]);
            _mm256_storeu_ps(&new_cont.label[i], _mm256_add_ps(v1, v2));
        }
#endif
    }
};

// DominanceFunction model

class dominance_spprc
{
public:
    inline bool operator()(const spp_spprc_res_cont& res_cont_1,
        const spp_spprc_res_cont& res_cont_2) const
    {
#ifndef _USE_SIMD_
        for (int i = 0; i < res_cont_1.label.size(); ++i) {
            if (res_cont_1.label[i] > res_cont_2.label[i]) {
                return false;
            }
        }
#else
        for (int i = 0; i < res_cont_1.label.size(); i += 8) {
            __m256 v1 = _mm256_loadu_ps(&res_cont_1.label[i]);
            __m256 v2 = _mm256_loadu_ps(&res_cont_2.label[i]);
            __m256 v3 = _mm256_cmp_ps(v1, v2, _CMP_GT_OS);
            if (!_mm256_testz_ps(v3, v3)) {
                return false;
            }
        }
#endif
        return true;
    }
};

#ifdef _REAL_DATA_
bool read_txt(std::vector<std::vector<float>>& result, std::string path) {
    if (!result.empty()) {
        return false;
    }
    std::ifstream file(path);
    if (file.is_open()) {
        int num_cols = 1;
        std::string line;
        std::getline(file, line); // skip first row
        std::istringstream iss(line);
        std::string value;
        while (iss >> value) {
            ++num_cols;
            if (iss.peek() == '\t')
                iss.ignore();
        }
        while (std::getline(file, line)) {
            std::vector<float> row;
            std::istringstream iss(line);
            float value;
            while (iss >> value) {
                row.push_back(value);

                if (iss.peek() == '\t')
                    iss.ignore();
                while (iss.peek() == '\t') {
                    row.push_back(-1);
                    iss.ignore();
                }
            }
            while (row.size() < num_cols) {
                row.push_back(-1);
            }
            result.push_back(row);
        }
        file.close();
        return true;
    }
    else {
        return false;
    }
}

int main()
{
    SPPRC_Example_Graph g;

    std::vector<std::vector<float>> net;
    std::vector<std::vector<float>> node;
    if (!read_txt(net, "net.txt")) {
        std::cout << "reading net.txt error" << std::endl;
        return -1;
    }
    if (!read_txt(node, "node.txt")) {
        std::cout << "reading node.txt error" << std::endl;
        return -1;
    }

    for (auto& row : node) {
        add_vertex(SPPRC_Example_Graph_Vert_Prop((int)row[0]), g);
    }

    int edge_i = 0;
    for (auto& row : net) {
#ifndef _USE_SIMD_
        std::vector<float> vec({ row[3], row[8], row[11], row[12], row[13] });
#else
        std::vector<float> vec({ row[3], row[8], row[11], row[12], row[13], (float)0, (float)0, (float)0 });
#endif
        add_edge((int)row[0], (int)row[1], SPPRC_Example_Graph_Arc_Prop(++edge_i, vec), g);
    }

    std::cout << "num of node:" << node.size() << std::endl;
    std::cout << "num of average out edge:" << net.size() * 1.0 / node.size() << std::endl;
    std::cout << "num of label:" << 5 << std::endl;

    boost::graph_traits< SPPRC_Example_Graph >::vertex_descriptor s = 1;
    boost::graph_traits< SPPRC_Example_Graph >::vertex_descriptor t = 12982;

    std::vector<
        std::vector< boost::graph_traits< SPPRC_Example_Graph >::edge_descriptor > >
        opt_solutions_spprc;
    std::vector< spp_spprc_res_cont > pareto_opt_rcs_spprc;
#ifndef _USE_SIMD_
    std::vector<float> init_vec(5, 0);
#else
    std::vector<float> init_vec(8, 0);
#endif
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "start..." << std::endl;
    r_c_shortest_paths(g, get(&SPPRC_Example_Graph_Vert_Prop::num, g),
        get(&SPPRC_Example_Graph_Arc_Prop::num, g), s, t, opt_solutions_spprc,
        pareto_opt_rcs_spprc, spp_spprc_res_cont(init_vec), ref_spprc(),
        dominance_spprc(),
        std::allocator< boost::r_c_shortest_paths_label< SPPRC_Example_Graph,
        spp_spprc_res_cont > >());
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    std::cout << "SPPRC:" << std::endl;
    std::cout << "Time: " << duration.count() / 1000000.0 << " second" << std::endl;
    std::cout << "Number of optimal solutions: ";
    std::cout << static_cast<int>(opt_solutions_spprc.size()) << std::endl;
    //for (int i = 0; i < static_cast<int>(opt_solutions_spprc.size()); ++i)
    //{
    //    std::cout << "Path: ";
    //    for (int j = static_cast<int>(opt_solutions_spprc[i].size()) - 1;
    //        j >= 0; --j)
    //        std::cout << source(opt_solutions_spprc[i][j], g) << " ";
    //    std::cout << num_node - 1 << std::endl;
    //}
    std::getchar();
    return 0;
}

#else
int main()
{
    SPPRC_Example_Graph g;

    int num_node = 2000;
    int avg_num_edge = 10;
    int num_label = 30;

    std::cout << "num of node:" << num_node << std::endl;
    std::cout << "num of average out edge:" << avg_num_edge << std::endl;
    std::cout << "num of label:" << num_label << std::endl;

    for (int i = 0; i < num_node; ++i) {
        add_vertex(SPPRC_Example_Graph_Vert_Prop(i), g);
    }

    std::default_random_engine e;
    std::uniform_int_distribution<int> u_num_out_node(1, 2 * avg_num_edge - 1);
    std::uniform_int_distribution<unsigned> u_out_node(0, num_node - 1);
    std::uniform_real_distribution<float> u_label(0, 100);
    int edge_i = 0;
    for (int i = 0; i < num_node; ++i) {
        for (int j = 0; j < u_num_out_node(e); ++j) {
            int tmp = u_out_node(e);
            while (tmp == i) {
                tmp = u_out_node(e);
            }
            std::vector<float> vec;
            for (int k = 0; k < num_label; ++k) {
                vec.push_back(u_label(e));
            }
#ifdef _USE_SIMD_
            // padding for AVX register
            for (int k = 0; k < (8 - num_label % 8) % 8; ++k) {
                vec.push_back((float)0);
            }
#endif
            add_edge(i, tmp, SPPRC_Example_Graph_Arc_Prop(i, vec), g);
        }
    }
    
    boost::graph_traits< SPPRC_Example_Graph >::vertex_descriptor s = 0;
    boost::graph_traits< SPPRC_Example_Graph >::vertex_descriptor t = num_node - 1;

    std::vector<
        std::vector< boost::graph_traits< SPPRC_Example_Graph >::edge_descriptor > >
        opt_solutions_spprc;
    std::vector< spp_spprc_res_cont > pareto_opt_rcs_spprc;
#ifndef _USE_SIMD_
    std::vector<float> init_vec(num_label, 0);
#else
    std::vector<float> init_vec(num_label + (8 - num_label % 8) % 8, 0);
#endif
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "start..." << std::endl;
    r_c_shortest_paths(g, get(&SPPRC_Example_Graph_Vert_Prop::num, g),
        get(&SPPRC_Example_Graph_Arc_Prop::num, g), s, t, opt_solutions_spprc,
        pareto_opt_rcs_spprc, spp_spprc_res_cont(init_vec), ref_spprc(),
        dominance_spprc(),
        std::allocator< boost::r_c_shortest_paths_label< SPPRC_Example_Graph,
        spp_spprc_res_cont > >());
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    std::cout << "SPPRC:" << std::endl;
    std::cout << "Time: " << duration.count() / 1000000.0 << " second" << std::endl;
    std::cout << "Number of optimal solutions: ";
    std::cout << static_cast<int>(opt_solutions_spprc.size()) << std::endl;
    //for (int i = 0; i < static_cast<int>(opt_solutions_spprc.size()); ++i)
    //{
    //    std::cout << "Path: ";
    //    for (int j = static_cast<int>(opt_solutions_spprc[i].size()) - 1;
    //        j >= 0; --j)
    //        std::cout << source(opt_solutions_spprc[i][j], g) << " ";
    //    std::cout << num_node - 1 << std::endl;
    //}
    std::getchar();
    return 0;
}
#endif // _REAL_DATA_





