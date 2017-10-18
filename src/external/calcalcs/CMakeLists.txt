add_library(pismcalcalcs OBJECT calcalcs.c utCalendar2_cal.c)
set_target_properties(pismcalcalcs PROPERTIES COMPILE_FLAGS "-w")
