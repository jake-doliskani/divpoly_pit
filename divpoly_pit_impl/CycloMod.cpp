#include "CycloMod.h"

CycloMod::CycloMod(long modulusDegree) {
    this->modulusDegree = modulusDegree;
}

CycloMod::~CycloMod() {
}

void CycloMod::mulMod(ZZ_pEX& result, 
                          const ZZ_pEX& f, 
                          const ZZ_pEX& g) {
    mul(result, f, g);
    reduce(result, result);
}

void CycloMod::sqrMod(ZZ_pEX& result, 
                          const ZZ_pEX& f) {
    sqr(result, f);
    reduce(result, result);
}

void CycloMod::reduce(ZZ_pEX& result, const ZZ_pEX& a) {
    if (deg(a) < modulusDegree) {
        result = a;
        return;
    }
    
    long degree = deg(a);
    long j = 0;
    ZZ_pEX temp;
    clear(temp);

    for (long i = 0; i <= degree; i++) {
        SetCoeff(temp, j, coeff(a, i) + coeff(temp, j));
        j++;
        
        if (j == modulusDegree)
            j = 0;
    }
    
    result = temp;
    temp.kill();
}

/**
 * Computes {@code f^e}. The code is basically a copy of NTL's powerMod.
 * @param result
 * @param f
 * @param e
 */
void CycloMod::powerMod(ZZ_pEX& result, const ZZ_pEX& f, const ZZ& e) {

   if (e == 0) {
      set(result);
      return;
   }

   if (e == 1) {
      result = f;
      return;
   }

   if (e == 2) {
      sqrMod(result, f);
      return;
   }

   long n = NumBits(e);

   ZZ_pEX res;
   res.SetMaxLength(modulusDegree);
   set(res);

   long i;

   if (n < 16) {
      // plain square-and-multiply algorithm

      for (i = n - 1; i >= 0; i--) {
         sqrMod(res, res);
         if (bit(e, i))
            mulMod(res, res, f);
      }

      result = res;
      return;
   }

   long k = optimalWinSize(n);
   k = min(k, 3);

   vec_ZZ_pEX v;

   v.SetLength(1L << (k-1));

   v[0] = f;
 
   if (k > 1) {
      ZZ_pEX t;
      sqrMod(t, f);

      for (i = 1; i < (1L << (k-1)); i++)
         mulMod(v[i], v[i-1], t);
   }


   long val;
   long cnt;
   long m;

   val = 0;
   for (i = n-1; i >= 0; i--) {
      val = (val << 1) | bit(e, i); 
      if (val == 0)
         sqrMod(res, res);
      else if (val >= (1L << (k-1)) || i == 0) {
         cnt = 0;
         while ((val & 1) == 0) {
            val = val >> 1;
            cnt++;
         }

         m = val;
         while (m > 0) {
            sqrMod(res, res);
            m = m >> 1;
         }

         mulMod(res, res, v[val >> 1]);

         while (cnt > 0) {
            sqrMod(res, res);
            cnt--;
         }

         val = 0;
      }
   }

   result = res;
}


long CycloMod::optimalWinSize(long n) {
   long k;
   double v, v_new;


   v = n/2.0 + 1.0;
   k = 1;

   for (;;) {
      v_new = n/(double(k+2)) + double(1L << k);
      if (v_new >= v) break;
      v = v_new;
      k++;
   }

   return k;
}