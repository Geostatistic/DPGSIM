/*
Implementação do plugin Distributed Independent Random Fields SGSIM.

(c) 2014, LPM/UFRGS, Luiz Gustavo Rasera, Péricles Lopes Machado

BASED ON AR2GEMS SGSIM IMPLEMENTATION
*/



#ifndef PLUGINS_LPM_UFRGS_COMMON_H_
#define PLUGINS_LPM_UFRGS_COMMON_H_

#include <string>

#if defined(_WIN32) || defined(WIN32)
  #ifdef LIB_STATIC
    #define PLUGINS_LPM_UFRGS_DECL
  #else
    #ifdef PLUGINS_LPM_UFRGS_EXPORT
      #define PLUGINS_LPM_UFRGS_DECL __declspec(dllexport)
    #else
      #define PLUGINS_LPM_UFRGS_DECL __declspec(dllimport)
    #endif
  #endif
#else
    #define PLUGINS_LPM_UFRGS_DECL
#endif

#include <cstring>

namespace LPM_UFRGS {
	inline std::string to_string(int i)
	{
		char buffer[100];
		sprintf(buffer, "%07d", i);
		
		return buffer;
	}
}



#endif // PLUGINS_LPM_UFRGS_COMMON_H_
