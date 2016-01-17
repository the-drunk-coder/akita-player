#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/program_options.hpp>
#include <iostream>

namespace iproc = boost::interprocess;
namespace po = boost::program_options;


// initialize the command line parameter parser
po::options_description init_opts(int ac, char *av[], po::variables_map *vm) {

  po::options_description desc("Message Queue Parameters");
  desc.add_options()
    ("help", "Display this help!")
    ("q", po::value<std::string>(), "A command queue!")
    ("msg", po::value<char>(), "A message to send!")
    ;
  // ----- end options ... what kind of syntax is this ??

  po::store(po::command_line_parser(ac, av).options(desc).run(),*vm);
  po::notify(*vm);

  return desc;
}

int main(int ac, char *av[])
{

  // initialize command line options
  po::variables_map vm;
  po::options_description desc = init_opts(ac, av, &vm);

  std::string cmd_queue_name = "";

  char msg;
  
  if (vm.count("q")) {
    cmd_queue_name = vm["q"].as<std::string>();
    if(cmd_queue_name.size() == 0){
      std::cout << "Please enter command queue name !" << std::endl;
      return 0;
    }
    
  }
  bool use_stdin = true;
  if (vm.count("msg")) {
    msg = vm["msg"].as<char>();
    use_stdin = false;  
  }


  
  try{
      
    //Create a message_queue.
    iproc:: message_queue out_q
      (iproc::open_or_create               //only create
       ,cmd_queue_name.c_str()           //name
       ,100                       //max message number
       ,sizeof(char)               //max message size
       );

    if(use_stdin){
      while(true){
	std::cin >> msg;
	out_q.send(&msg, sizeof(msg), 0);
	 
      } 
    } else {
      out_q.send(&msg, sizeof(msg), 0);
    }
      
  }
  catch(iproc::interprocess_exception &ex){
    std::cout << ex.what() << std::endl;
    return 1;
  }

  return 0;
}
