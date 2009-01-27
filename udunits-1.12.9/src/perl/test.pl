#!/usr/local/bin/perl


use UDUNITS;

$epsilon = 1e-5;

print STDERR "Initializing units module......................";
UDUNITS::init("../lib/udunits.dat") == 0 || die "error\n";
print STDERR "ok\n";

print STDERR "Getting blank unit.............................";
$unit = UDUNITS::new();
print STDERR "ok\n";

print STDERR "Verifying blank unit...........................";
!$unit->hasorigin() || die "error\n";
!$unit->istime() || die "error\n";
$spec = $unit->print();
$spec =~ /^$/ || die "error\n";
print STDERR "ok\n";

print STDERR "Scanning ratio units...........................";
$megaparsec_barn = UDUNITS::scan("megaparsec barn") || die "error\n";
$tsp = UDUNITS::scan("tsp") || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying ratio unit...........................";
!$tsp->hasorigin() || die "error\n";
!$tsp->istime() || die "error\n";
print STDERR "ok\n";

print STDERR "Converting ratio units.........................";
$megaparsec_barn->convert($tsp, $slope, $intercept) == 0
    || die "Couldn't convert units\n";
abs(1 - $slope/0.626035) < $epsilon && $intercept == 0 || 
    die "error\n";
print STDERR "ok\n";

print STDERR "Scanning temporal unit.........................";
$moment = UDUNITS::scan("seconds since 1925-01-01 00:00:00") || 
    die "error\n";
print STDERR "ok\n";

print STDERR "Verifying temporal unit........................";
$moment->istime() || die "error\n";
print STDERR "ok\n";

print STDERR "Printing temporal unit.........................";
$timestring = $moment->print() || die "error\n";
$timestring =~ /seconds since 1925-01-01 00:00:00.0* UTC/ || die "error\n";
print STDERR "ok\n";

print STDERR "Converting value to calendar time..............";
$moment->valtocal(5, $year, $month, $day, $hour, $minute, $second) == 0 || 
    die "error\n";
$year == 1925 && $month == 1 && $day == 1 && $hour == 0 && $minute == 0
    && abs(1 - $second/5) < $epsilon || die "error\n";
print STDERR "ok\n";

print STDERR "Converting calendar time to value..............";
$offset = $moment->caltoval($year, $month, $day, $hour, $minute, $second);
abs(1 - $offset/5) < $epsilon || die "error\n";
print STDERR "ok\n";

print STDERR "Scanning interval units........................";
$celsius = UDUNITS::scan("Celsius @ 100") || die "error\n";
$fahrenheit = UDUNITS::scan("Fahrenheit") || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying interval unit........................";
$celsius->hasorigin() || die "error\n";
print STDERR "ok\n";

print STDERR "Converting interval units......................";
$celsius->convert($fahrenheit, $slope, $intercept) == 0
    || die "Couldn't convert units\n";
abs(1 - $slope/1.8) < $epsilon && abs(1 - $intercept/212) < $epsilon ||
    die "error\n";
print STDERR "ok\n";

print STDERR "Verifying clear()..............................";
$unit = UDUNITS::scan("second");
$unit->clear();
$spec = $unit->print();
$spec =~ /^$/ || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying shift()..............................";
$unit = UDUNITS::scan("meter");
$unit->shift(2);
$spec = $unit->print() || die "error\n";
$spec =~ /meter @ 2/ || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying scale()..............................";
$unit = UDUNITS::scan("meter");
$unit->scale(2);
$spec = $unit->print() || die "error\n";
#print STDERR "spec = $spec\n";
$spec =~ /2 meter/ || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying multiply()...........................";
$unit1 = UDUNITS::scan("meter");
$unit2 = UDUNITS::scan("second");
$unit1->multiply($unit2);
$spec = $unit1->print() || die "error\n";
#print STDERR "spec = $spec\n";
$spec =~ /meter second/ || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying invert().............................";
$unit = UDUNITS::scan("meter");
$unit->invert();
$spec = $unit->print() || die "error\n";
#print STDERR "spec = $spec\n";
$spec =~ /meter-1/ || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying divide().............................";
$unit1 = UDUNITS::scan("meter");
$unit2 = UDUNITS::scan("second");
$unit1->divide($unit2);
$spec = $unit1->print() || die "error\n";
#print STDERR "spec = $spec\n";
$spec =~ /meter second-1/ || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying raise()..............................";
$unit = UDUNITS::scan("meter");
$unit->raise(-2);
$spec = $unit->print() || die "error\n";
#print STDERR "spec = $spec\n";
$spec =~ /meter-2/ || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying dup()................................";
$unit1 = UDUNITS::scan("meter");
$unit2 = $unit1->dup();
$unit1->convert($unit2, $slope, $intercept);
abs(1 - $slope) < $epsilon && $intercept == 0 || die "error\n";
print STDERR "ok\n";

print STDERR "Verifying UDUNITS::dup().......................";
$unit1 = UDUNITS::scan("meter");
$unit2 = utUnitPtr::dup($unit1);
$unit1->convert($unit2, $slope, $intercept);
abs(1 - $slope) < $epsilon && $intercept == 0 || die "error\n";
print STDERR "ok\n";

print STDERR "Terminating units module.......................";
UDUNITS::term();
print STDERR "ok\n";
