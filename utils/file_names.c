#include <dirent.h> 
#include <stdlib.h>
#include <string.h>

void* open_dir(const char * const dir_name){
  DIR* d;
  d = opendir(".");
  return d;
}

void close_dir(DIR *d){
  if (d) closedir(d);
}

void next_file(DIR *d, char ch[256], int* len){
  do {
    struct dirent *dir = readdir(d);
    if (dir){
      if (dir->d_type == DT_REG)
      {
        strncpy(ch, dir->d_name, 256);
        size_t slen = strlen(dir->d_name);
        if (slen > 255) {
          *len = 255;
        } else { 
          *len = (int) slen;          
        }
      } else {
        *len = -1;
      }
    } else {
      *len = 0;
    }
  }
  while (*len<0);  
}
