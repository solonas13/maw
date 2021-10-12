#ifndef __STREAM_H_INCLUDED
#define __STREAM_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <string>
#include <bitset>

#include "./utils.h"
#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

template<typename T>
struct stream_reader {
  stream_reader(std::string fname, INT buf_bytes = (4L << 20))
      : m_bufelems(buf_bytes / sizeof(T)),
        m_filled(0),
        m_offset(0L) {
    m_file = utils::file_open(fname, "r");
	m_buffer = (T*) malloc(m_bufelems* sizeof(T));
    refill();
  }

  ~stream_reader() {
    free( m_buffer);
    std::fclose(m_file);
  }

  inline T peek() {
    return m_buffer[m_pos];
  }

  inline T& operator[] (INT i) {
    assert(i >= m_offset);
    if (i >= m_offset + m_filled) {  // elem not in a buffer
      if (i < m_offset + m_filled + m_bufelems) {
        // elem in the next buffer, no need to fseek
        refill();
      } else {
        // elem potentially far away, fseek
        refill(i);
      }
    }

    return m_buffer[i - m_offset];
  }

  inline T getValue(INT i) {
    assert(i >= m_offset);
    if (i >= m_offset + m_filled) {  // elem not in a buffer
      if (i < m_offset + m_filled + m_bufelems) {
        // elem in the next buffer, no need to fseek
        refill();
      } else {
        // elem potentially far away, fseek
        refill(i);
      }
    }
    if (i<m_offset)
	refill(i);

    return m_buffer[i - m_offset];
  }

  stream_reader& operator++ () {
    ++m_pos;
    if (m_pos == m_filled) refill();

    return *this;
  }

  inline bool goto_set(){
    std::fseek(m_file, 0, SEEK_SET);
    m_offset=0;
    m_filled=0;
    refill();
    return true;
  }

  inline bool goto_pos(INT offset){
    m_filled=0;
    refill(offset);
    return true;
  }

  inline bool goto_end(INT n_elems){
    std::fseek(m_file, 0, SEEK_END);
    m_offset=n_elems;
    m_filled=0;
    refill_reverse(1);
    return true;
  }

  inline T read() {
    T ret = m_buffer[m_pos++];
    if (m_pos == m_filled) refill();

    return ret;
  }

  inline bool empty() {
    return (!m_filled && !refill());
  }

  inline INT getpos(){
    return this->m_pos;
   }
  inline void setpos(INT pos){
    m_pos=pos; 
   }

   inline T read_reverse(bool B=false){
     if (B) refill_reverse(1);
     T ret = m_buffer[m_pos--];
     if (m_pos==-1) refill_reverse(2);
     return ret;
   }

  inline bool empty_reverse() {
    return (!m_filled && !refill_reverse(2));
  }



 private:
  INT refill(INT new_offset = -1) {
    if (new_offset != -1) {
      std::fseek(m_file, new_offset, SEEK_SET);
      m_offset=new_offset/sizeof(T);
      if (new_offset != m_offset * sizeof(T))
	std::cout<<" il y a un pti pb ici"<<std::endl;
    } else {
      m_offset += m_filled;
    }

    m_filled = (INT)std::fread(m_buffer, sizeof(T), m_bufelems, m_file);
    m_pos = 0;

    return m_filled;
  }

  INT refill_reverse(INT time) {
    bool B=true;
    if (ftell(m_file)>= m_bufelems*sizeof(T)*time)
    	std::fseek(m_file, -m_bufelems*sizeof(T)*time, SEEK_CUR);
    else
    {
        std::fseek(m_file, 0, SEEK_SET);
        B=false;
    }
    
    if (B)
        m_filled = (INT)std::fread(m_buffer, sizeof(T), m_bufelems, m_file);
    else
        m_filled = (INT)std::fread(m_buffer, sizeof(T), m_offset, m_file);

    if (m_offset==0)
    {
        m_filled=0;
	return 0;
    }
    m_offset -= m_filled;
    m_pos = m_filled-1;

    return m_filled;
  }


  std::FILE *m_file;

  INT m_bufelems, m_filled, m_pos;
  INT m_offset;
  T *m_buffer;
};

template<typename T>
struct stream_writer {
  stream_writer(std::string fname, INT bufsize = (4 << 20), int a=0)
      : m_bufelems(bufsize / sizeof(T)),
        m_filled(0) {
    if (a==0)
    {	m_file = utils::file_open(fname.c_str(), "w");}
    else
    {
	if ( ! ( m_file = fopen ( fname.c_str(), "a") ) )
        {
                fprintf ( stderr, " Error: Cannot open file %s!\n", fname.c_str() );
        }
     }
     m_buffer =(T*) malloc(m_bufelems*sizeof(T));

  }

  ~stream_writer() {
    if (m_filled) flush();
    free( m_buffer);
    std::fclose(m_file);
  }

  inline void write(T x) {
    m_buffer[m_filled++] = x;
    if (m_filled == m_bufelems) flush();
  }

  inline INT getpos(){
	flush();
	return ftell(m_file);
 }

 private:
  void flush() {
    utils::add_objects_to_file(m_file, m_filled, m_buffer);
    m_filled = 0;
  }

  std::FILE *m_file;
  INT m_bufelems, m_filled;
  T *m_buffer;
};

 
#endif  // __STREAM_H_INCLUDED
