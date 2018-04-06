add_library(calcalcs OBJECT calcalcs.c utCalendar2_cal.c)
set_target_properties(calcalcs PROPERTIES COMPILE_FLAGS "-w")
