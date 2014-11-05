#include "log.hpp"

void initLog() {
  logging::add_file_log (
    keywords::file_name = "myLog.log",                                        /*< file name pattern >*/
    keywords::format = (
      // This makes the sink to write log records that look like this:
      // YYYY-MM-DD HH:MI:SS: <normal> A normal severity message
      expr::stream
      << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
      << ": <" << logging::trivial::severity
      << "> " << expr::smessage
      )
    );

  logging::core::get()->set_filter (
    logging::trivial::severity >= logging::trivial::info
    );
}