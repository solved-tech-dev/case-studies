// Reject ordinary file writes during the compute interval, including the
// initial GPU setup. Device handles and MPI shared memory are not output files.
#include <atomic>
#include <cstdarg>
#include <cstring>
#include <dlfcn.h>
#include <fcntl.h>
#include <cerrno>
#include <cstdio>
static std::atomic<bool> active{false};
static std::atomic<unsigned long> rejected{0};
static char firstPath[512]={};
extern "C" void benchIoStart(){rejected=0;firstPath[0]=0;active=true;}
extern "C" unsigned long benchIoStop(){active=false;return rejected.load();}
extern "C" const char* benchIoFirstPath(){return firstPath;}
static bool deny(const char* p,int f){
 if(!active||!(f&(O_WRONLY|O_RDWR|O_CREAT|O_TRUNC|O_APPEND)))return false;
 if(p&&(std::strncmp(p,"/dev/",5)==0||std::strncmp(p,"/proc/",6)==0||std::strncmp(p,"/sys/",5)==0))return false;
 if(rejected.fetch_add(1)==0&&p)std::strncpy(firstPath,p,sizeof(firstPath)-1);
 errno=EPERM;return true;
}
extern "C" int open(const char* p,int f,...){
 mode_t mode=0;if(f&O_CREAT){va_list ap;va_start(ap,f);mode=va_arg(ap,int);va_end(ap);}
 if(deny(p,f))return -1;
 static auto fn=reinterpret_cast<int(*)(const char*,int,...)>(dlsym(RTLD_NEXT,"open"));return fn(p,f,mode);
}
extern "C" int open64(const char* p,int f,...){
 mode_t mode=0;if(f&O_CREAT){va_list ap;va_start(ap,f);mode=va_arg(ap,int);va_end(ap);}
 if(deny(p,f))return -1;
 static auto fn=reinterpret_cast<int(*)(const char*,int,...)>(dlsym(RTLD_NEXT,"open64"));return fn(p,f,mode);
}
extern "C" int openat(int fd,const char* p,int f,...){
 mode_t mode=0;if(f&O_CREAT){va_list ap;va_start(ap,f);mode=va_arg(ap,int);va_end(ap);}
 if(deny(p,f))return -1;
 static auto fn=reinterpret_cast<int(*)(int,const char*,int,...)>(dlsym(RTLD_NEXT,"openat"));return fn(fd,p,f,mode);
}
extern "C" FILE* fopen(const char* p,const char* mode){
 if(deny(p,(std::strpbrk(mode,"wa+")?O_WRONLY:O_RDONLY)))return nullptr;
 static auto fn=reinterpret_cast<FILE*(*)(const char*,const char*)>(dlsym(RTLD_NEXT,"fopen"));return fn(p,mode);
}
extern "C" FILE* fopen64(const char* p,const char* mode){
 if(deny(p,(std::strpbrk(mode,"wa+")?O_WRONLY:O_RDONLY)))return nullptr;
 static auto fn=reinterpret_cast<FILE*(*)(const char*,const char*)>(dlsym(RTLD_NEXT,"fopen64"));return fn(p,mode);
}
