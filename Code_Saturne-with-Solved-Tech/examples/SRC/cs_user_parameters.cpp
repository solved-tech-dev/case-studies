#include "cs_headers.h"
#include <cstdlib>
#include <dlfcn.h>

extern "C" void cs_user_linear_solvers(void) {
  const char *library = std::getenv("CS_LINGLONGAMG_ADAPTER_LIBRARY");
  if (library == nullptr)
    bft_error(__FILE__, __LINE__, 0, "CS_LINGLONGAMG_ADAPTER_LIBRARY is required");
  // Keep the module loaded until SLES finalization invokes its callbacks.
  void *module = dlopen(library, RTLD_NOW | RTLD_GLOBAL);
  if (module == nullptr)
    bft_error(__FILE__, __LINE__, 0, "Cannot load LinglongAMG adapter: %s", dlerror());
  auto define = reinterpret_cast<void (*)(int)>(dlsym(module, "cs_linglongamg_define"));
  if (define == nullptr)
    bft_error(__FILE__, __LINE__, 0, "Missing cs_linglongamg_define");
  define(CS_F_(p)->id);
}
