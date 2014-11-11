#include "log.hpp"

void initLog() {

  // adds LineID and TimeStamp attributes
  logging::add_common_attributes();

  // define formatter for logs. 'auto' keyword avoids compile error
  // explanation at http://stackoverflow.com/questions/17766998/boost-log-why-doesnt-this-compile
  auto fmt = expr::stream
    << expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
    << ": <" << logging::trivial::severity << "> " << expr::smessage;

  logging::add_file_log(keywords::file_name = "dynamix.log",
    keywords::auto_flush = true, keywords::format = fmt);

  logging::add_console_log(std::clog,
    keywords::auto_flush = true, keywords::format = fmt);

  logging::core::get()->set_filter(
    logging::trivial::severity >= logging::trivial::info);

  return;
}