#!/bin/sh

cd mr
autoreconf -i
./configure --prefix=$PWD/build

# patch MR
if [ `md5 -q mr/tsil.cpp` = 'f8ffff19291c1561b473a463cd5c994e' ]; then
    patch -p1 <<_END_
--- a/mr/tsil.cpp
+++ b/mr/tsil.cpp
@@ -34,12 +34,12 @@ namespace mr

   std::complex<long double> Li2(std::complex<long double> z)
   {
-    return TSIL_Dilog(z.real()+1.0I*z.imag());
+    return c2pp(TSIL_Dilog(z.real()+1.0I*z.imag()));
   }

   std::complex<long double> Li3(std::complex<long double> z)
   {
-    return TSIL_Trilog(z.real()+1.0I*z.imag());
+    return c2pp(TSIL_Trilog(z.real()+1.0I*z.imag()));
   }

   std::complex<long double> acc(std::complex<long double> z)
_END_
fi

make install -j6
