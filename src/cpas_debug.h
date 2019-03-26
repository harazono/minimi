#ifndef _DEBUG_HEADER
#define _DEBUG_HEADER

#include <cstring>

#define DEBUGMACRO_FUNCTIONNAME __PRETTY_FUNCTION__
#define DEBUGMACRO_ABORT() GDB_On_SEGV::gdb_handler(0)
#ifndef DEBUGMACRO_ERRORSTREAM
 #define DEBUGMACRO_ERRORSTREAM std::cerr
#endif

inline void DebugOutputVariableDumpHelper1(){}
inline void DebugOutputVariableDumpHelper1(const char* str){if(std::strlen(str) != 0){DEBUGMACRO_ERRORSTREAM << (str) << " = ";}}
inline void DebugOutputVariableDumpHelper3(){}
inline void DebugOutputVariableDumpHelper3(const char* str){if(std::strlen(str) != 0){DEBUGMACRO_ERRORSTREAM << "              " << (str) << " = ";}}
inline void DebugOutputVariableDumpHelper2(){}
template<class T> inline void DebugOutputVariableDumpHelper2(const T& arg){DEBUGMACRO_ERRORSTREAM << (arg) << std::endl;}

#define VARDUMP(arg, ...)   (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),VARDUMP_1(__VA_ARGS__))
#define VARDUMP_1(arg, ...) (DebugOutputVariableDumpHelper3(#arg),DebugOutputVariableDumpHelper2(arg),VARDUMP_2(__VA_ARGS__))
#define VARDUMP_2(arg, ...) (DebugOutputVariableDumpHelper3(#arg),DebugOutputVariableDumpHelper2(arg),VARDUMP_3(__VA_ARGS__))
#define VARDUMP_3(arg, ...) (DebugOutputVariableDumpHelper3(#arg),DebugOutputVariableDumpHelper2(arg),VARDUMP_4(__VA_ARGS__))
#define VARDUMP_4(arg, ...) (DebugOutputVariableDumpHelper3(#arg),DebugOutputVariableDumpHelper2(arg),VARDUMP_5(__VA_ARGS__))
#define VARDUMP_5(arg, ...) (DebugOutputVariableDumpHelper3(#arg),DebugOutputVariableDumpHelper2(arg),VARDUMP_6(__VA_ARGS__))
#define VARDUMP_6(arg, ...) (DebugOutputVariableDumpHelper3(#arg),DebugOutputVariableDumpHelper2(arg),VARDUMP_7(__VA_ARGS__))
#define VARDUMP_7(arg, ...) (DebugOutputVariableDumpHelper3(#arg),DebugOutputVariableDumpHelper2(arg),VARDUMP_8(__VA_ARGS__))
#define VARDUMP_8(arg, ...) (DebugOutputVariableDumpHelper3(#arg),DebugOutputVariableDumpHelper2(arg),VARDUMP_9(__VA_ARGS__))
#define VARDUMP_9(arg, ...) (DebugOutputVariableDumpHelper3(#arg),DebugOutputVariableDumpHelper2(arg))
#define DUMP(arg, ...)   (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),DUMP_1(__VA_ARGS__))
#define DUMP_1(arg, ...) (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),DUMP_2(__VA_ARGS__))
#define DUMP_2(arg, ...) (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),DUMP_3(__VA_ARGS__))
#define DUMP_3(arg, ...) (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),DUMP_4(__VA_ARGS__))
#define DUMP_4(arg, ...) (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),DUMP_5(__VA_ARGS__))
#define DUMP_5(arg, ...) (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),DUMP_6(__VA_ARGS__))
#define DUMP_6(arg, ...) (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),DUMP_7(__VA_ARGS__))
#define DUMP_7(arg, ...) (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),DUMP_8(__VA_ARGS__))
#define DUMP_8(arg, ...) (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg),DUMP_9(__VA_ARGS__))
#define DUMP_9(arg, ...) (DebugOutputVariableDumpHelper1(#arg),DebugOutputVariableDumpHelper2(arg))

#if defined(NDEBUG)
#define MYASSERT_WMD(message, condition, dumpmethod) ;
#else
#define MYASSERT_WMD(message, condition, dumpmethod) if(condition) ; else {\
  DEBUGMACRO_ERRORSTREAM << "ASSERTION FAILED!\n" \
  "  File      : " << __FILE__ << "\n" \
  "  Line      : " << __LINE__ << "\n" \
  "  Function  : " << DEBUGMACRO_FUNCTIONNAME << "\n" \
  "  Condition : " << #condition << "\n" \
  "  Message   : " << message << "\n" \
  "  Dump      : "; \
  dumpmethod; DEBUGMACRO_ABORT(); \
  }
#endif

#define MYASSERT_WM(message, condition)                    MYASSERT_WMD(message, condition, DEBUGMACRO_ERRORSTREAM << std::flush)
#define MYASSERT(condition)                        MYASSERT_WM("N/A", condition)

#define MYASSERT_EQUALS_WMD(message, arg1, arg2, dumpmethod) MYASSERT_WMD(message, (arg1) == (arg2), (VARDUMP(arg1, arg2), (dumpmethod)))
#define MYASSERT_EQUALS_WM(message, arg1, arg2)              MYASSERT_EQUALS_WMD(message, arg1, arg2, DEBUGMACRO_ERRORSTREAM << std::flush)
#define MYASSERT_EQUALS(arg1, arg2)                          MYASSERT_EQUALS_WM("N/A", arg1, arg2)

#define MYASSERT_LT_WMD(message, arg1, arg2, dumpmethod)     MYASSERT_WMD(message, (arg1) < (arg2), (VARDUMP(arg1, arg2), (dumpmethod)))
#define MYASSERT_LT_WM(message, arg1, arg2)                  MYASSERT_LT_WMD(message, arg1, arg2, DEBUGMACRO_ERRORSTREAM << std::flush)
#define MYASSERT_LT(arg1, arg2)                              MYASSERT_LT_WM("N/A", arg1, arg2)

#define MYASSERT_LTE_WMD(message, arg1, arg2, dumpmethod)    MYASSERT_WMD(message, (arg1) <= (arg2), (VARDUMP(arg1, arg2), (dumpmethod)))
#define MYASSERT_LTE_WM(message, arg1, arg2)                 MYASSERT_LTE_WMD(message, arg1, arg2, DEBUGMACRO_ERRORSTREAM << std::flush)
#define MYASSERT_LTE(arg1, arg2)                             MYASSERT_LTE_WM("N/A", arg1, arg2)

#define MYASSERT_GT_WMD(message, arg1, arg2, dumpmethod)     MYASSERT_WMD(message, (arg1) > (arg2), (VARDUMP(arg1, arg2), (dumpmethod)))
#define MYASSERT_GT_WM(message, arg1, arg2)                  MYASSERT_GT_WMD(message, arg1, arg2, DEBUGMACRO_ERRORSTREAM << std::flush)
#define MYASSERT_GT(arg1, arg2)                              MYASSERT_GT_WM("N/A", arg1, arg2)

#define MYASSERT_GTE_WMD(message, arg1, arg2, dumpmethod)    MYASSERT_WMD(message, (arg1) >= (arg2), (VARDUMP(arg1, arg2), (dumpmethod)))
#define MYASSERT_GTE_WM(message, arg1, arg2)                 MYASSERT_GTE_WMD(message, arg1, arg2, DEBUGMACRO_ERRORSTREAM << std::flush)
#define MYASSERT_GTE(arg1, arg2)                             MYASSERT_GTE_WM("N/A", arg1, arg2)

#define MYASSERT_NEVERREACH_WD(dumpmethod)        MYASSERT_WMD("Logic error", false, dumpmethod)
#define MYASSERT_NEVERREACH( )            MYASSERT_NEVERREACH_WD(DEBUGMACRO_ERRORSTREAM << std::flush)

#define MYASSERT_NOTIMPLEMENTED_WD(dumpmethod)    MYASSERT_WMD("Not implemented", false, dumpmethod)
#define MYASSERT_NOTIMPLEMENTED( ) MYASSERT_NEVERREACH_WD(DEBUGMACRO_ERRORSTREAM << std::flush)

#endif // #ifndef _DEBUG_HEADER

