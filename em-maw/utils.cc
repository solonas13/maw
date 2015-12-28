#include "./utils.h"

#include <sys/time.h>
#include <errno.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>
#include <fstream>


namespace utils {

long double wclock() {
  timeval tim;
  gettimeofday(&tim, NULL);

  return tim.tv_sec + (tim.tv_usec / 1000000.0L);
}

void file_copy(std::string src_fname, std::string dest_fname) {
  std::ifstream src(src_fname.c_str(), std::ios::binary);
  std::ofstream dst(dest_fname.c_str(), std::ios::binary);
  dst << src.rdbuf();
}

void file_delete(std::string fname) {
  int res = std::remove(fname.c_str());
  if (res) {
    fprintf(stderr, "Failed to delete %s: %s\n",
        fname.c_str(), strerror(errno));
    std::exit(EXIT_FAILURE);
  }
}

std::string absolute_path(std::string fname) {
  char path[1 << 18];
  bool created = false;

  if (!file_exists(fname)) {
    // We need to create the file, since realpath fails on non-existing files.
    std::fclose(file_open(fname, "w"));
    created = true;
  }
  if (!realpath(fname.c_str(), path)) {
    fprintf(stderr, "Error: realpath failed for %s\n", fname.c_str());
    std::exit(EXIT_FAILURE);
  }

  if (created)
    file_delete(fname);

  return std::string(path);
}

std::FILE *file_open(std::string fname, std::string mode) {
  std::FILE *f = std::fopen(fname.c_str(), mode.c_str());
  if (!f) {
    std::perror(fname.c_str());
    std::exit(EXIT_FAILURE);
  }

  return f;
}

long file_size(std::string fname) {
  std::FILE *f = file_open(fname, "rt");
  std::fseek(f, 0L, SEEK_END);
  long size = std::ftell(f);
  std::fclose(f);

  return size;
}

bool file_exists(std::string fname) {
  std::FILE *f = std::fopen(fname.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL)
    std::fclose(f);

  return ret;
}

void find_stxxl_config() {
  if (file_exists("./.stxxl")) {
    fprintf(stderr, "STXXL config file detected.\n");
    return;
  } else if (file_exists(std::string(std::getenv("HOME")) + "/.stxxl")) {
    fprintf(stderr, "Cannot find STXXL config file. Using $HOME/.stxxl\n");
    std::string src = std::string(std::getenv("HOME")) + "/.stxxl";
    std::string dest = absolute_path("./.stxxl");
    file_copy(src, dest);
    return;
  } else {
    fprintf(stderr, "Error: failed to find/copy STXXL config file!\n");
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace utils
