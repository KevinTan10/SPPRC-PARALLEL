// modify from r_c_shortest_paths.hpp

#ifndef R_C_SHORTEST_PATHS_PARALLEL_HPP
#define R_C_SHORTEST_PATHS_PARALLEL_HPP

#include <map>
#include <queue>
#include <vector>
#include <list>

#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/property_map/property_map.hpp>

#include <mutex>
#include <condition_variable>
#include <thread>
#include <iostream>
#include <atomic>


namespace boost
{
    // r_c_shortest_paths_label struct
    template < class Graph, class Resource_Container >
    struct r_c_shortest_paths_label
        : public boost::enable_shared_from_this<
        r_c_shortest_paths_label< Graph, Resource_Container > >
    {
        r_c_shortest_paths_label(const unsigned long n,
            const Resource_Container& rc = Resource_Container(),
            const boost::shared_ptr<
            r_c_shortest_paths_label< Graph, Resource_Container > >
            pl
            = boost::shared_ptr<
            r_c_shortest_paths_label< Graph, Resource_Container > >(),
            const typename graph_traits< Graph >::edge_descriptor& ed
            = graph_traits< Graph >::edge_descriptor(),
            const typename graph_traits< Graph >::vertex_descriptor& vd
            = graph_traits< Graph >::vertex_descriptor())
            : num(n)
            , cumulated_resource_consumption(rc)
            , p_pred_label(pl)
            , pred_edge(ed)
            , resident_vertex(vd)
            , b_is_dominated(false)
        {
        }

        r_c_shortest_paths_label& operator=(const r_c_shortest_paths_label& other)
        {
            if (this == &other)
                return *this;
            this->~r_c_shortest_paths_label();
            new (this) r_c_shortest_paths_label(other);
            return *this;
        }
        const unsigned long num;
        Resource_Container cumulated_resource_consumption;
        const boost::shared_ptr<
            r_c_shortest_paths_label< Graph, Resource_Container > >
            p_pred_label;
        const typename graph_traits< Graph >::edge_descriptor pred_edge;
        const typename graph_traits< Graph >::vertex_descriptor resident_vertex;
        bool b_is_dominated;
    }; // r_c_shortest_paths_label

    template < class Graph, class Resource_Container >
    inline bool operator==(
        const r_c_shortest_paths_label< Graph, Resource_Container >& l1,
        const r_c_shortest_paths_label< Graph, Resource_Container >& l2)
    {
        return l1.cumulated_resource_consumption
            == l2.cumulated_resource_consumption;
    }

    template < class Graph, class Resource_Container >
    inline bool operator!=(
        const r_c_shortest_paths_label< Graph, Resource_Container >& l1,
        const r_c_shortest_paths_label< Graph, Resource_Container >& l2)
    {
        return !(l1 == l2);
    }

    template < class Graph, class Resource_Container >
    inline bool operator<(
        const r_c_shortest_paths_label< Graph, Resource_Container >& l1,
        const r_c_shortest_paths_label< Graph, Resource_Container >& l2)
    {
        return l1.cumulated_resource_consumption
            < l2.cumulated_resource_consumption;
    }

    template < class Graph, class Resource_Container >
    inline bool operator>(
        const r_c_shortest_paths_label< Graph, Resource_Container >& l1,
        const r_c_shortest_paths_label< Graph, Resource_Container >& l2)
    {
        return l2.cumulated_resource_consumption
            < l1.cumulated_resource_consumption;
    }

    template < class Graph, class Resource_Container >
    inline bool operator<=(
        const r_c_shortest_paths_label< Graph, Resource_Container >& l1,
        const r_c_shortest_paths_label< Graph, Resource_Container >& l2)
    {
        return l1 < l2 || l1 == l2;
    }

    template < class Graph, class Resource_Container >
    inline bool operator>=(
        const r_c_shortest_paths_label< Graph, Resource_Container >& l1,
        const r_c_shortest_paths_label< Graph, Resource_Container >& l2)
    {
        return l2 < l1 || l1 == l2;
    }

    template < typename Graph, typename Resource_Container >
    inline bool operator<(
        const boost::shared_ptr<
        r_c_shortest_paths_label< Graph, Resource_Container > >& t,
        const boost::shared_ptr<
        r_c_shortest_paths_label< Graph, Resource_Container > >& u)
    {
        return *t < *u;
    }

    template < typename Graph, typename Resource_Container >
    inline bool operator<=(
        const boost::shared_ptr<
        r_c_shortest_paths_label< Graph, Resource_Container > >& t,
        const boost::shared_ptr<
        r_c_shortest_paths_label< Graph, Resource_Container > >& u)
    {
        return *t <= *u;
    }

    template < typename Graph, typename Resource_Container >
    inline bool operator>(
        const boost::shared_ptr<
        r_c_shortest_paths_label< Graph, Resource_Container > >& t,
        const boost::shared_ptr<
        r_c_shortest_paths_label< Graph, Resource_Container > >& u)
    {
        return *t > *u;
    }

    template < typename Graph, typename Resource_Container >
    inline bool operator>=(
        const boost::shared_ptr<
        r_c_shortest_paths_label< Graph, Resource_Container > >& t,
        const boost::shared_ptr<
        r_c_shortest_paths_label< Graph, Resource_Container > >& u)
    {
        return *t >= *u;
    }

    namespace detail
    {

        // r_c_shortest_paths_dispatch function (body/implementation)
        template < class Graph, class VertexIndexMap, class EdgeIndexMap,
            class Resource_Container, class Resource_Extension_Function,
            class Dominance_Function, class Label_Allocator >
        void r_c_shortest_paths_dispatch(const Graph& g,
            const VertexIndexMap& vertex_index_map,
            const EdgeIndexMap& /*edge_index_map*/,
            typename graph_traits< Graph >::vertex_descriptor s,
            typename graph_traits< Graph >::vertex_descriptor t,
            // each inner vector corresponds to a pareto-optimal path
            std::vector<
            std::vector< typename graph_traits< Graph >::edge_descriptor > >&
            pareto_optimal_solutions,
            std::vector< Resource_Container >& pareto_optimal_resource_containers,
            bool b_all_pareto_optimal_solutions,
            // to initialize the first label/resource container
            // and to carry the type information
            const Resource_Container& rc, Resource_Extension_Function& ref,
            Dominance_Function& dominance,
            // to specify the memory management strategy for the labels
            Label_Allocator /*la*/)
        {
            pareto_optimal_resource_containers.clear();
            pareto_optimal_solutions.clear();

            size_t i_label_num = 0;
#if defined(BOOST_NO_CXX11_ALLOCATOR)
            typedef typename Label_Allocator::template rebind<
                r_c_shortest_paths_label< Graph, Resource_Container > >::other
                LAlloc;
#else
            typedef typename std::allocator_traits< Label_Allocator >::
                template rebind_alloc<
                r_c_shortest_paths_label< Graph, Resource_Container > >
                LAlloc;
            typedef std::allocator_traits< LAlloc > LTraits;
#endif
            LAlloc l_alloc;
            typedef boost::shared_ptr<
                r_c_shortest_paths_label< Graph, Resource_Container > >
                Splabel;
            std::queue< Splabel >
                unprocessed_labels; // queue instead of priority queue for multiple label

            bool b_feasible = true;
            Splabel splabel_first_label = boost::allocate_shared<
                r_c_shortest_paths_label< Graph, Resource_Container > >(l_alloc,
                    i_label_num++, rc,
                    boost::shared_ptr<
                    r_c_shortest_paths_label< Graph, Resource_Container > >(),
                    typename graph_traits< Graph >::edge_descriptor(), s);

            unprocessed_labels.push(splabel_first_label);
            std::vector< std::list< Splabel > > vec_vertex_labels_data(
                num_vertices(g));
            iterator_property_map<
                typename std::vector< std::list< Splabel > >::iterator,
                VertexIndexMap >
                vec_vertex_labels(vec_vertex_labels_data.begin(), vertex_index_map);
            vec_vertex_labels[s].push_back(splabel_first_label);
            typedef std::vector< typename std::list< Splabel >::iterator >
                vec_last_valid_positions_for_dominance_data_type;
            vec_last_valid_positions_for_dominance_data_type
                vec_last_valid_positions_for_dominance_data(num_vertices(g));
            iterator_property_map<
                typename vec_last_valid_positions_for_dominance_data_type::iterator,
                VertexIndexMap >
                vec_last_valid_positions_for_dominance(
                    vec_last_valid_positions_for_dominance_data.begin(),
                    vertex_index_map);
            BGL_FORALL_VERTICES_T(v, g, Graph)
            {
                put(vec_last_valid_positions_for_dominance, v,
                    vec_vertex_labels[v].begin());
            }
            std::vector< size_t > vec_last_valid_index_for_dominance_data(
                num_vertices(g), 0);
            iterator_property_map< std::vector< size_t >::iterator, VertexIndexMap >
                vec_last_valid_index_for_dominance(
                    vec_last_valid_index_for_dominance_data.begin(),
                    vertex_index_map);
            std::vector< bool > b_vec_vertex_already_checked_for_dominance_data(
                num_vertices(g), false);
            iterator_property_map< std::vector< bool >::iterator, VertexIndexMap >
                b_vec_vertex_already_checked_for_dominance(
                    b_vec_vertex_already_checked_for_dominance_data.begin(),
                    vertex_index_map);

            // vertex label buffer
            std::vector< std::list< Splabel > > vec_vertex_labels_buffer_data(
                num_vertices(g));
            iterator_property_map<
                typename std::vector< std::list< Splabel > >::iterator,
                VertexIndexMap >
                vec_vertex_labels_buffer(vec_vertex_labels_buffer_data.begin(), vertex_index_map);

            // lock
            std::vector< std::mutex > vec_vertex_labels_mutex_data(
                num_vertices(g));
            iterator_property_map< std::vector< std::mutex >::iterator, VertexIndexMap >
                vec_vertex_labels_mutex(vec_vertex_labels_mutex_data.begin(), vertex_index_map);
            std::vector< std::mutex > vec_vertex_labels_buffer_mutex_data(
                num_vertices(g));
            iterator_property_map< std::vector< std::mutex >::iterator, VertexIndexMap >
                vec_vertex_labels_buffer_mutex(vec_vertex_labels_buffer_mutex_data.begin(), vertex_index_map);
            std::mutex queue_mutex;
            std::condition_variable condition;

            // multithread
            bool stop = 0;
            std::atomic<size_t> active_thread_num(0);
            std::vector<std::thread> threads;
            for (size_t i = 0; i < 32; ++i) {
                threads.emplace_back([&] {
                    std::unique_lock<std::mutex> queue_lock(queue_mutex, std::defer_lock);
                    while (true) {
                        queue_lock.lock();
                        condition.wait(queue_lock, [&] { return stop || !unprocessed_labels.empty(); });

                        if (stop) return;

                        ++active_thread_num;
                        Splabel cur_label = unprocessed_labels.front();
                        unprocessed_labels.pop();
                        queue_lock.unlock();

                        std::unique_lock<std::mutex> vertex_lock(get(vec_vertex_labels_mutex, cur_label->resident_vertex));

                        // flush buffer
                        {
                            std::unique_lock<std::mutex> vertex_buffer_lock(get(vec_vertex_labels_buffer_mutex, cur_label->resident_vertex));
                            std::list< Splabel >& list_labels_cur_vertex
                                = get(vec_vertex_labels, cur_label->resident_vertex);
                            list_labels_cur_vertex.splice(list_labels_cur_vertex.end(), get(vec_vertex_labels_buffer, cur_label->resident_vertex));
                        }

                        if (!cur_label->b_is_dominated)
                        {
                            typename boost::graph_traits< Graph >::vertex_descriptor
                                i_cur_resident_vertex
                                = cur_label->resident_vertex;
                            std::list< Splabel >& list_labels_cur_vertex
                                = get(vec_vertex_labels, i_cur_resident_vertex);
                            if (vec_last_valid_index_for_dominance[i_cur_resident_vertex]
                                < list_labels_cur_vertex.size()
                                && list_labels_cur_vertex.size() >= 2)
                            {
                                typename std::list< Splabel >::iterator outer_iter
                                    = list_labels_cur_vertex.begin();
                                bool b_outer_iter_at_or_beyond_last_valid_pos_for_dominance
                                    = false;
                                while (outer_iter != list_labels_cur_vertex.end())
                                {
                                    Splabel& cur_outer_splabel = *outer_iter;
                                    typename std::list< Splabel >::iterator inner_iter
                                        = outer_iter;
                                    if (!b_outer_iter_at_or_beyond_last_valid_pos_for_dominance
                                        && outer_iter
                                        == get(vec_last_valid_positions_for_dominance,
                                            i_cur_resident_vertex))
                                        b_outer_iter_at_or_beyond_last_valid_pos_for_dominance
                                        = true;
                                    if (!get(b_vec_vertex_already_checked_for_dominance,
                                        i_cur_resident_vertex)
                                        || b_outer_iter_at_or_beyond_last_valid_pos_for_dominance)
                                    {
                                        ++inner_iter;
                                    }
                                    else
                                    {
                                        inner_iter
                                            = get(vec_last_valid_positions_for_dominance,
                                                i_cur_resident_vertex);
                                        ++inner_iter;
                                    }
                                    bool b_outer_iter_erased = false;
                                    while (inner_iter != list_labels_cur_vertex.end())
                                    {
                                        Splabel& cur_inner_splabel = *inner_iter;
                                        if (dominance(cur_outer_splabel
                                            ->cumulated_resource_consumption,
                                            cur_inner_splabel
                                            ->cumulated_resource_consumption))
                                        {
                                            cur_inner_splabel->b_is_dominated = true;
                                            cur_inner_splabel.reset();
                                            typename std::list< Splabel >::iterator buf
                                                = inner_iter;
                                            ++inner_iter;
                                            list_labels_cur_vertex.erase(buf);
                                            continue;
                                        }
                                        else
                                            ++inner_iter;
                                        if (dominance(cur_inner_splabel
                                            ->cumulated_resource_consumption,
                                            cur_outer_splabel
                                            ->cumulated_resource_consumption))
                                        {
                                            cur_outer_splabel->b_is_dominated = true;
                                            cur_outer_splabel.reset();
                                            typename std::list< Splabel >::iterator buf
                                                = outer_iter;
                                            ++outer_iter;
                                            list_labels_cur_vertex.erase(buf);
                                            b_outer_iter_erased = true;
                                            break;
                                        }
                                    }
                                    if (!b_outer_iter_erased)
                                        ++outer_iter;
                                }
                                if (list_labels_cur_vertex.size() > 1)
                                    put(vec_last_valid_positions_for_dominance,
                                        i_cur_resident_vertex,
                                        (--(list_labels_cur_vertex.end())));
                                else
                                    put(vec_last_valid_positions_for_dominance,
                                        i_cur_resident_vertex,
                                        list_labels_cur_vertex.begin());
                                put(b_vec_vertex_already_checked_for_dominance,
                                    i_cur_resident_vertex, true);
                                put(vec_last_valid_index_for_dominance,
                                    i_cur_resident_vertex,
                                    list_labels_cur_vertex.size() - 1);
                            }
                        }
                        if (!cur_label->b_is_dominated)
                        {
                            typename graph_traits< Graph >::vertex_descriptor cur_vertex
                                = cur_label->resident_vertex;
                            typename graph_traits< Graph >::out_edge_iterator oei, oei_end;
                            std::vector< Splabel > new_labels;
                            for (boost::tie(oei, oei_end) = out_edges(cur_vertex, g);
                                oei != oei_end; ++oei)
                            {
                                Splabel new_label = boost::allocate_shared<
                                    r_c_shortest_paths_label< Graph, Resource_Container > >(
                                        l_alloc, i_label_num++,
                                        cur_label->cumulated_resource_consumption, cur_label,
                                        *oei, target(*oei, g));
                                ref(g,
                                    new_label->cumulated_resource_consumption,
                                    new_label->p_pred_label->cumulated_resource_consumption,
                                    new_label->pred_edge);
                                {
                                    std::unique_lock<std::mutex> vertex_buffer_lock(get(vec_vertex_labels_buffer_mutex, new_label->resident_vertex));
                                    vec_vertex_labels_buffer[new_label->resident_vertex].push_back(
                                        new_label);
                                }
                                new_labels.push_back(new_label);
                            }
                            queue_lock.lock();
                            for (auto& new_label : new_labels) {
                                unprocessed_labels.push(new_label);
                            }
                            queue_lock.unlock();
                            condition.notify_all();
                        }
                        else
                        {
                            cur_label.reset();
                        }
                        --active_thread_num;
                    }
                });
            }

            while (active_thread_num != 0 || !unprocessed_labels.empty()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            stop = true;
            condition.notify_all();
            for (size_t i = 0; i < threads.size(); ++i) {
                threads[i].join();
            }

            std::list< Splabel > dsplabels = get(vec_vertex_labels, t);
            typename std::list< Splabel >::const_iterator csi = dsplabels.begin();
            typename std::list< Splabel >::const_iterator csi_end = dsplabels.end();
            // if d could be reached from o
            if (!dsplabels.empty())
            {
                for (; csi != csi_end; ++csi)
                {
                    std::vector< typename graph_traits< Graph >::edge_descriptor >
                        cur_pareto_optimal_path;
                    boost::shared_ptr<
                        r_c_shortest_paths_label< Graph, Resource_Container > >
                        p_cur_label = *csi;
                    pareto_optimal_resource_containers.push_back(
                        p_cur_label->cumulated_resource_consumption);
                    while (p_cur_label->num != 0)
                    {
                        cur_pareto_optimal_path.push_back(p_cur_label->pred_edge);
                        p_cur_label = p_cur_label->p_pred_label;
                    }
                    pareto_optimal_solutions.push_back(cur_pareto_optimal_path);
                    if (!b_all_pareto_optimal_solutions)
                        break;
                }
            }

            BGL_FORALL_VERTICES_T(i, g, Graph)
            {
                std::list< Splabel >& list_labels_cur_vertex = vec_vertex_labels[i];
                typename std::list< Splabel >::iterator si
                    = list_labels_cur_vertex.begin();
                const typename std::list< Splabel >::iterator si_end
                    = list_labels_cur_vertex.end();
                for (; si != si_end; ++si)
                {
                    (*si).reset();
                }
            }
        } // r_c_shortest_paths_dispatch

    } // detail


    // r_c_shortest_paths functions (handle/interface)
    // - return all pareto-optimal solutions
    // - specify Label_Allocator
    template < class Graph, class VertexIndexMap, class EdgeIndexMap,
        class Resource_Container, class Resource_Extension_Function,
        class Dominance_Function, class Label_Allocator >
    void r_c_shortest_paths(const Graph& g, const VertexIndexMap& vertex_index_map,
        const EdgeIndexMap& edge_index_map,
        typename graph_traits< Graph >::vertex_descriptor s,
        typename graph_traits< Graph >::vertex_descriptor t,
        // each inner vector corresponds to a pareto-optimal path
        std::vector<
        std::vector< typename graph_traits< Graph >::edge_descriptor > >&
        pareto_optimal_solutions,
        std::vector< Resource_Container >& pareto_optimal_resource_containers,
        // to initialize the first label/resource container
        // and to carry the type information
        const Resource_Container& rc, const Resource_Extension_Function& ref,
        const Dominance_Function& dominance,
        // to specify the memory management strategy for the labels
        Label_Allocator la)
    {
        r_c_shortest_paths_dispatch(g, vertex_index_map, edge_index_map, s, t,
            pareto_optimal_solutions, pareto_optimal_resource_containers, true, rc,
            ref, dominance, la);
    }
    // r_c_shortest_paths

    // check_r_c_path function
    template < class Graph, class Resource_Container,
        class Resource_Extension_Function >
    void check_r_c_path(const Graph& g,
        const std::vector< typename graph_traits< Graph >::edge_descriptor >&
        ed_vec_path,
        const Resource_Container& initial_resource_levels,
        // if true, computed accumulated final resource levels must
        // be equal to desired_final_resource_levels
        // if false, computed accumulated final resource levels must
        // be less than or equal to desired_final_resource_levels
        bool b_result_must_be_equal_to_desired_final_resource_levels,
        const Resource_Container& desired_final_resource_levels,
        Resource_Container& actual_final_resource_levels,
        const Resource_Extension_Function& ref, bool& b_is_a_path_at_all,
        bool& b_feasible, bool& b_correctly_extended,
        typename graph_traits< Graph >::edge_descriptor& ed_last_extended_arc)
    {
        size_t i_size_ed_vec_path = ed_vec_path.size();
        std::vector< typename graph_traits< Graph >::edge_descriptor > buf_path;
        if (i_size_ed_vec_path == 0)
            b_feasible = true;
        else
        {
            if (i_size_ed_vec_path == 1
                || target(ed_vec_path[0], g) == source(ed_vec_path[1], g))
                buf_path = ed_vec_path;
            else
                for (size_t i = i_size_ed_vec_path; i > 0; --i)
                    buf_path.push_back(ed_vec_path[i - 1]);
            for (size_t i = 0; i < i_size_ed_vec_path - 1; ++i)
            {
                if (target(buf_path[i], g) != source(buf_path[i + 1], g))
                {
                    b_is_a_path_at_all = false;
                    b_feasible = false;
                    b_correctly_extended = false;
                    return;
                }
            }
        }
        b_is_a_path_at_all = true;
        b_feasible = true;
        b_correctly_extended = false;
        Resource_Container current_resource_levels = initial_resource_levels;
        actual_final_resource_levels = current_resource_levels;
        for (size_t i = 0; i < i_size_ed_vec_path; ++i)
        {
            ed_last_extended_arc = buf_path[i];
            b_feasible = ref(g, actual_final_resource_levels,
                current_resource_levels, buf_path[i]);
            current_resource_levels = actual_final_resource_levels;
            if (!b_feasible)
                return;
        }
        if (b_result_must_be_equal_to_desired_final_resource_levels)
            b_correctly_extended
            = actual_final_resource_levels == desired_final_resource_levels
            ? true
            : false;
        else
        {
            if (actual_final_resource_levels < desired_final_resource_levels
                || actual_final_resource_levels == desired_final_resource_levels)
                b_correctly_extended = true;
        }
    } // check_path

} // namespace

#endif // BOOST_GRAPH_R_C_SHORTEST_PATHS_HPP