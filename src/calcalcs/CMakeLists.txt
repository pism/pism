add_library(pismcalcalcs calcalcs.c utCalendar2_cal.c)
target_link_libraries (pismcalcalcs ${UDUNITS2_LIBRARY})
set_target_properties(pismcalcalcs PROPERTIES COMPILE_FLAGS "-w")

