#include "./include/mapreduce.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sys/time.h>

namespace wordcount {

class map_task;
class reduce_task;

typedef
mapreduce::job<
  wordcount::map_task,
  wordcount::reduce_task>
job;

class map_task : boost::noncopyable
{
  public:
    typedef std::string   key_type;
    typedef std::ifstream value_type;
    typedef std::string   intermediate_key_type;
    typedef unsigned      intermediate_value_type;

    map_task(job::map_task_runner &runner)
      : runner_(runner)
    {
    }

    // 'value_type' is not a reference to const to enable streams to be passed
    //    key: input filename
    //    value: ifstream of the open file
    void operator()(key_type const &/*key*/, value_type &value)
    {
        while (!value.eof())
        {
            std::string word;
            value >> word;
            std::transform(word.begin(), word.end(), word.begin(),
                           std::bind1st(
                               std::mem_fun(&std::ctype<char>::tolower),
                               &std::use_facet<std::ctype<char> >(std::locale::classic())));

            runner_.emit_intermediate(word, 1);
        }
    }

  private:
    job::map_task_runner &runner_;
};

class reduce_task : boost::noncopyable
{
  public:
    typedef std::string  key_type;
    typedef size_t       value_type;

    reduce_task(job::reduce_task_runner &runner)
      : runner_(runner)
    {
    }

    template<typename It>
    void operator()(typename map_task::intermediate_key_type const &key, It it, It ite)
    {
        reduce_task::value_type result = 0;
        for (; it!=ite; ++it)
           result += *it;
        runner_.emit(key, result);
    }

  private:
    job::reduce_task_runner &runner_;
};

}   // namespace wordcount

int main(int argc, char **argv)
{
    wordcount::job::datasource_type datasource;
    datasource.set_directory(argv[1]);

    mapreduce::specification  spec;
    mapreduce::results        result;
    wordcount::job            mr(datasource);
    mr.run<mapreduce::schedule_policy::cpu_parallel>(spec, result);

    // output the results
    std::cout << std::endl << "\n" << "MapReduce statistics:";
    std::cout << "\n  " << "MapReduce job runtime                     : " << result.job_runtime << " seconds, of which...";
    std::cout << "\n  " << "  Map phase runtime                       : " << result.map_runtime << " seconds";
    std::cout << "\n  " << "  Reduce phase runtime                    : " << result.reduce_runtime << " seconds";
    std::cout << "\n\n  " << "Map:";
    std::cout << "\n    " << "Total Map keys                          : " << result.counters.map_tasks;
    std::cout << "\n    " << "Map keys processed                      : " << result.counters.map_tasks_completed;
    std::cout << "\n    " << "Map key processing errors               : " << result.counters.map_tasks_error;
    std::cout << "\n    " << "Number of Map Tasks run (in parallel)   : " << result.counters.actual_map_tasks;
    std::cout << "\n    " << "Fastest Map key processed in            : " << *std::min_element(result.map_times.begin(), result.map_times.end()) << " seconds";
    std::cout << "\n    " << "Slowest Map key processed in            : " << *std::max_element(result.map_times.begin(), result.map_times.end()) << " seconds";
    std::cout << "\n    " << "Average time to process Map keys        : " << std::accumulate(result.map_times.begin(), result.map_times.end(), boost::int64_t()) / result.map_times.size() << " seconds";
    std::cout << "\n\n  " << "Reduce:";
    std::cout << "\n    " << "Number of Reduce Tasks run (in parallel): " << result.counters.actual_reduce_tasks;
    std::cout << "\n    " << "Number of Result Files                  : " << result.counters.num_result_files;
    std::cout << "\n    " << "Fastest Reduce key processed in         : " << *std::min_element(result.reduce_times.begin(), result.reduce_times.end()) << " seconds";
    std::cout << "\n    " << "Slowest Reduce key processed in         : " << *std::max_element(result.reduce_times.begin(), result.reduce_times.end()) << " seconds";
    std::cout << "\n    " << "Average time to process Reduce keys     : " << std::accumulate(result.reduce_times.begin(), result.reduce_times.end(), boost::int64_t()) / result.map_times.size() << " seconds";

    return 0;
}