#!/usr/bin/awk -f

BEGIN{c=1}
{
if ($13!="H") {
  printf( \
    "%-6s  %3d  %-4s%3s   %3d     %7.3f %7.3f %7.3f  %4.2f  %4.2f     %4s%4s\n",\
    $1,c,$3,$4,$6,$7,$8,$9+z,$10,$11,$12,$13);
    c+=1
  }
}