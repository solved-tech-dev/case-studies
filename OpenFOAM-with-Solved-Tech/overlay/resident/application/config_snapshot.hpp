// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cstdlib>
#include <cstring>
#include <climits>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdexcept>
#include <initializer_list>

inline bool resident_configuration_snapshots() {
 const char* enabled=std::getenv("RESIDENT_CONFIG_PRELOADED");
 if(!enabled)return false;
 if(std::strcmp(enabled,"1"))throw std::runtime_error("invalid configuration preload mode");
 const char prefix[]="/proc/self/fd/";
 for(const char* name:{"RESIDENT_GS_CONFIG","RESIDENT_PRESSURE_CONFIG","RESIDENT_PROVIDER_CONFIG"}) {
  const char* path=std::getenv(name);
  if(!path||std::strncmp(path,prefix,sizeof(prefix)-1))throw std::runtime_error("preloaded configuration descriptor required");
  char* end=nullptr;const long value=std::strtol(path+sizeof(prefix)-1,&end,10);
  if(end==path+sizeof(prefix)-1||*end||value<0||value>INT_MAX)throw std::runtime_error("invalid configuration descriptor");
  const int fd=static_cast<int>(value);struct stat info{};char target[256];
  const auto length=::readlink(path,target,sizeof(target)-1);
  if(::fstat(fd,&info)||!S_ISREG(info.st_mode)||info.st_size<=0||length<=0||length>=static_cast<ssize_t>(sizeof(target)-1))throw std::runtime_error("invalid configuration snapshot");
  target[length]=0;
  const int seals=F_SEAL_WRITE|F_SEAL_GROW|F_SEAL_SHRINK|F_SEAL_SEAL;
  const char memory_prefix[]="/memfd:xk2647-";const int actual=fcntl(fd,F_GET_SEALS);
  if(std::strncmp(target,memory_prefix,sizeof(memory_prefix)-1)||actual<0||(actual&seals)!=seals)throw std::runtime_error("configuration snapshot must be sealed memory");
 }
 return true;
}

inline bool resident_diagnostic_enabled() {
 for(const char* name:{"RESIDENT_STATE_DUMPS","XK_PRIVATE_PRESSURE_CAPTURE","RESIDENT_CUPTI_TRANSFER_LOG",
  "XK_PRIVATE_RESIDENT_PROFILE","XK_SDK_PROFILE","XK_PRIVATE_RESIDENT_SYSTEM_FINGERPRINT",
  "XK_PRIVATE_RESIDENT_POINTER_TRACE","XK_PRIVATE_RESIDENT_STORAGE_LOG","XK_PRIVATE_RESIDENT_MEMORY_LOG"}) {
  const char* setting=std::getenv(name);if(setting&&*setting)return true;
 }
 return false;
}
