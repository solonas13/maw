#ifndef __UTILS_H_INCLUDED
#define __UTILS_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
#include <sstream>


namespace utils {

#define STRX(x) #x
#define STR(x) STRX(x)

long double wclock();

std::string absolute_path(std::string fname);
std::FILE *file_open(std::string fname, std::string mode);
void file_copy(std::string src_fname, std::string dest_fname);
void file_delete(std::string fname);
long file_size(std::string fname);
bool file_exists(std::string fname);
void find_stxxl_config();

template<typename T>
void add_objects_to_file(std::FILE *f, long count, T *tab) {
  long fwrite_ret = std::fwrite(tab, sizeof(T), count, f);
  if (fwrite_ret != count) {
    fprintf(stderr, "Error: fwrite in line %s of %s returned %ld\n",
        STR(__LINE__), STR(__FILE__), fwrite_ret);
    std::exit(EXIT_FAILURE);
  }
}

template<typename T>
void read_objects_from_file(std::FILE *f, long count, T *dest) {
  long fread_ret = std::fread(dest, sizeof(T), count, f);
  if (fread_ret != count) {
    fprintf(stderr, "Error: fread in line %s of %s returned %ld\n",
        STR(__LINE__), STR(__FILE__), fread_ret);
    std::exit(EXIT_FAILURE);
  }
}

template<typename T>
void read_objects_at_offset(std::FILE *f, long offset, long count, T *dest) {
  std::fseek(f, sizeof(T) * offset, SEEK_SET);
  read_objects_from_file<T>(f, count, dest);
}

template<typename T>
void read_objects_at_offset(std::string filename, long offset, long count,
    T *dest) {
  std::FILE *f = file_open(filename, "r");
  read_objects_at_offset<T>(f, offset, count, dest);
  std::fclose(f);
}

template<typename T>
T read_object_at_offset(std::FILE *f, long offset) {
  T ret;
  read_objects_at_offset(f, offset, 1L, &ret);

  return ret;
}

}  // namespace utils

#endif  // __UTILS_H_INCLUDED
