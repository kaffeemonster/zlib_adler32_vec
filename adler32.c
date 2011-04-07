/* adler32.c -- compute the Adler-32 checksum of a data stream
 * Copyright (C) 1995-2007 Mark Adler
 * For conditions of distribution and use, see copyright notice in zlib.h
 */

/* @(#) $Id$ */

#include "zutil.h"

#define local static

#define GCC_VERSION_GE(x) ((__GNUC__-0) * 100 + __GNUC_MINOR__-0 >= x)

#if GCC_VERSION_GE(301)
/* sometimes leakes out of old kernel header */
#  undef noinline
#  define noinline __attribute__((__noinline__))
#else
#  ifndef noinline
#    define noinline
#  endif
#endif

#if GCC_VERSION_GE(301)
# define GCC_ATTR_UNUSED_PARAM __attribute__((__unused__))
#else
# define GCC_ATTR_UNUSED_PARAM
#endif

#if GCC_VERSION_GE(296)
#  define likely(x)   __builtin_expect(!!(x), 1)
#  define unlikely(x) __builtin_expect(!!(x), 0)
#else
#  define likely(x)   (x)
#  define unlikely(x) (x)
#endif

#define ROUND_TO(x , n) ((x) & ~((n) - 1L))
#define DIV_ROUNDUP(a, b) (((a) + (b) - 1) / (b))
#define ALIGN_DIFF(x, n) (((intptr_t)((x)+(n) - 1L) & ~((intptr_t)(n) - 1L)) - (intptr_t)(x))
#define ALIGN_DOWN(x, n) (((intptr_t)(x)) & ~((intptr_t)(n) - 1L))
#define ALIGN_DOWN_DIFF(x, n) (((intptr_t)(x)) & ((intptr_t)(n) - 1L))

local uLong adler32_combine_(uLong adler1, uLong adler2, z_off64_t len2);

#define BASE 65521UL    /* largest prime smaller than 65536 */
#define NMAX 5552
/* NMAX is the largest n such that 255n(n+1)/2 + (n+1)(BASE-1) <= 2^32-1 */

#define DO1(buf,i)  {adler += (buf)[i]; sum2 += adler;}
#define DO2(buf,i)  DO1(buf,i); DO1(buf,i+1);
#define DO4(buf,i)  DO2(buf,i); DO2(buf,i+2);
#define DO8(buf,i)  DO4(buf,i); DO4(buf,i+4);
#define DO16(buf)   DO8(buf,0); DO8(buf,8);

#if defined(__alpha__)
/* even if gcc can generate a mul by inverse, the code is really
 * ugly (find global const pool pointer, load constant, a mul, lots
 * of shifts/add/sub), up to 14 instructions. The replacement code
 * only needs >= 5 instructions
 */
#  define NO_DIVIDE
#endif

/* use NO_DIVIDE if your processor does not do division in hardware */
#ifdef NO_DIVIDE
/* use NO_SHIFT if your processor does shift > 1 by loop */
#  ifdef NO_SHIFT
#    define reduce_full(a) \
    do { \
        if (a >= (BASE << 16)) a -= (BASE << 16); \
        if (a >= (BASE << 15)) a -= (BASE << 15); \
        if (a >= (BASE << 14)) a -= (BASE << 14); \
        if (a >= (BASE << 13)) a -= (BASE << 13); \
        if (a >= (BASE << 12)) a -= (BASE << 12); \
        if (a >= (BASE << 11)) a -= (BASE << 11); \
        if (a >= (BASE << 10)) a -= (BASE << 10); \
        if (a >= (BASE << 9)) a -= (BASE << 9); \
        if (a >= (BASE << 8)) a -= (BASE << 8); \
        if (a >= (BASE << 7)) a -= (BASE << 7); \
        if (a >= (BASE << 6)) a -= (BASE << 6); \
        if (a >= (BASE << 5)) a -= (BASE << 5); \
        if (a >= (BASE << 4)) a -= (BASE << 4); \
        if (a >= (BASE << 3)) a -= (BASE << 3); \
        if (a >= (BASE << 2)) a -= (BASE << 2); \
        if (a >= (BASE << 1)) a -= (BASE << 1); \
        if (a >= BASE) a -= BASE; \
    } while (0)
#    define reduce_x(a) \
    do { \
        if (MIN_WORK >= (1 << 6) && a >= (BASE << 6)) a -= (BASE << 6); \
        if (MIN_WORK >= (1 << 5) && a >= (BASE << 5)) a -= (BASE << 5); \
        if (a >= (BASE << 4)) a -= (BASE << 4); \
        if (a >= (BASE << 3)) a -= (BASE << 3); \
        if (a >= (BASE << 2)) a -= (BASE << 2); \
        if (a >= (BASE << 1)) a -= (BASE << 1); \
        if (a >= BASE) a -= BASE; \
    } while (0)
#    define reduce(a) reduce_full(a)
#  else
#    define reduce_full(a) \
    do { \
        unsigned long b = a & 0x0000ffff; \
        a >>= 16; \
        b -= a; \
        a <<= 4; \
        a += b; \
    } while(a >= BASE)
#    define reduce_x(a) \
    do { \
        unsigned long b = a & 0x0000ffff; \
        a >>= 16; \
        b -= a; \
        a <<= 4; \
        a += b; \
        a = a >= BASE ? a - BASE : a; \
    } while(0)
#    define reduce(a) \
    do { \
        unsigned long b = a & 0x0000ffff; \
        a >>= 16; \
        b -= a; \
        a <<= 4; \
        a += b; \
    } while(0)
#  endif
#else
#  define reduce_full(a) a %= BASE
#  define reduce_x(a) a %= BASE
#  define reduce(a) a %= BASE
#endif

local int host_is_bigendian()
{
    local const union {
        uInt d;
        unsigned char endian[sizeof(uInt)];
    } x = {1};
    return x.endian[0] == 0;
}

#ifndef NO_ADLER32_VEC
#  if defined(__arm__)
#    include "adler32_arm.c"
#  elif defined(__powerpc__) || defined(__powerpc64__)
#    include "adler32_ppc.c"
#  elif defined(__i386__) || defined(__x86_64__)
#    include "adler32_x86.c"
#  endif
#endif

#ifndef MIN_WORK
#  define MIN_WORK 16
#endif

/* ========================================================================= */
local noinline uLong adler32_1(adler, buf, len)
    uLong adler;
    const Bytef *buf;
    uInt len GCC_ATTR_UNUSED_PARAM;
{
    unsigned long sum2;

    /* split Adler-32 into component sums */
    sum2 = (adler >> 16) & 0xffff;
    adler &= 0xffff;

    adler += buf[0];
    if (adler >= BASE)
        adler -= BASE;
    sum2 += adler;
    if (sum2 >= BASE)
        sum2 -= BASE;
    return adler | (sum2 << 16);
}

/* ========================================================================= */
local noinline uLong adler32_common(adler, buf, len)
    uLong adler;
    const Bytef *buf;
    uInt len;
{
    unsigned long sum2;

    /* split Adler-32 into component sums */
    sum2 = (adler >> 16) & 0xffff;
    adler &= 0xffff;

    while (len--) {
        adler += *buf++;
        sum2 += adler;
    }
    if (adler >= BASE)
        adler -= BASE;
    reduce_x(sum2);             /* only added so many BASE's */
    return adler | (sum2 << 16);
}

#ifndef HAVE_ADLER32_VEC
#  if (defined(__LP64__) || ((SIZE_MAX-0) >> 31) >= 2) && !defined(NO_ADLER32_VEC)

/* On 64 Bit archs, we can do pseudo SIMD with a nice win.
 * This is esp. important for old Alphas, they do not have byte
 * access.
 * This needs some register but x86_64 is fine (>= 9 for the mainloop
 * req.). If your 64 Bit arch is more limited, throw it away...
 */
#    ifndef UINT64_C
#      if defined(_MSC_VER) || defined(__BORLANDC__)
#        define UINT64_C(c)    (c ## ui64)
#      else
#        define UINT64_C(c)    (c ## ULL)
#      endif
#    endif

#    undef VNMAX
#    define VNMAX (2*NMAX+((9*NMAX)/10))

/* ========================================================================= */
local noinline uLong adler32_vec(adler, buf, len)
    uLong adler;
    const Bytef *buf;
    uInt len;
{
    unsigned int s1, s2;
    unsigned int k;

    /* split Adler-32 into component sums */
    s1 = adler & 0xffff;
    s2 = (adler >> 16) & 0xffff;

    /* align input data */
    k    = ALIGN_DIFF(buf, sizeof(size_t));
    len -= k;
    if (k) do {
        s1 += *buf++;
        s2 += s1;
    } while(--k);

    k = len > VNMAX ? VNMAX : len;
    len -= k;
    if (likely(k >= 2 * sizeof(size_t))) do
    {
        unsigned int vs1, vs2;
        unsigned int vs1s;

        /* add s1 to s2 for rounds to come */
        s2 += s1 * ROUND_TO(k, sizeof(size_t));
        vs1s = vs1 = vs2 = 0;
        do {
            size_t vs1l = 0, vs1h = 0, vs1l_s = 0, vs1h_s = 0;
            unsigned int a, b, c, d, e, f, g, h;
            unsigned int j;

            j = k > 23 * sizeof(size_t) ? 23 : k/sizeof(size_t);
            k -= j * sizeof(size_t);
            /* add s1 to s1 round sum for rounds to come */
            vs1s += j * vs1;
            do {
                size_t in8 = *(const size_t *)buf;
                buf += sizeof(size_t);
                /* add this s1 to s1 round sum */
                vs1l_s += vs1l;
                vs1h_s += vs1h;
                /* add up input data to s1 */
                vs1l +=  in8 & UINT64_C(0x00ff00ff00ff00ff);
                vs1h += (in8 & UINT64_C(0xff00ff00ff00ff00)) >> 8;
            } while(--j);

            /* split s1 */
            if(host_is_bigendian()) {
                a = (vs1h >> 48) & 0x0000ffff;
                b = (vs1l >> 48) & 0x0000ffff;
                c = (vs1h >> 32) & 0x0000ffff;
                d = (vs1l >> 32) & 0x0000ffff;
                e = (vs1h >> 16) & 0x0000ffff;
                f = (vs1l >> 16) & 0x0000ffff;
                g = (vs1h      ) & 0x0000ffff;
                h = (vs1l      ) & 0x0000ffff;
            } else {
                a = (vs1l      ) & 0x0000ffff;
                b = (vs1h      ) & 0x0000ffff;
                c = (vs1l >> 16) & 0x0000ffff;
                d = (vs1h >> 16) & 0x0000ffff;
                e = (vs1l >> 32) & 0x0000ffff;
                f = (vs1h >> 32) & 0x0000ffff;
                g = (vs1l >> 48) & 0x0000ffff;
                h = (vs1h >> 48) & 0x0000ffff;
            }

            /* add s1 & s2 horiz. */
            vs2 += 8*a + 7*b + 6*c + 5*d + 4*e + 3*f + 2*g + 1*h;
            vs1 += a + b + c + d + e + f + g + h;

            /* split and add up s1 round sum */
            vs1l_s = ((vs1l_s      ) & UINT64_C(0x0000ffff0000ffff)) +
                     ((vs1l_s >> 16) & UINT64_C(0x0000ffff0000ffff));
            vs1h_s = ((vs1h_s      ) & UINT64_C(0x0000ffff0000ffff)) +
                     ((vs1h_s >> 16) & UINT64_C(0x0000ffff0000ffff));
            vs1l_s += vs1h_s;
            vs1s += ((vs1l_s      ) & UINT64_C(0x00000000ffffffff)) +
                    ((vs1l_s >> 32) & UINT64_C(0x00000000ffffffff));
        } while (k >= sizeof(size_t));
        reduce(vs1s);
        s2 += vs1s * 8 + vs2;
        reduce(s2);
        s1 += vs1;
        reduce(s1);
        len += k;
        k = len > VNMAX ? VNMAX : len;
        len -= k;
    } while (k >= sizeof(size_t));

    /* handle trailer */
    if (k) do {
        s1 += *buf++;
        s2 += s1;
    } while (--k);
    reduce(s1);
    reduce(s2);

    /* return recombined sums */
    return (s2 << 16) | s1;
}

#  else

/* ========================================================================= */
local noinline uLong adler32_vec(adler, buf, len)
    uLong adler;
    const Bytef *buf;
    uInt len;
{
    unsigned long sum2;
    unsigned n;

    /* split Adler-32 into component sums */
    sum2 = (adler >> 16) & 0xffff;
    adler &= 0xffff;

    /* do length NMAX blocks -- requires just one modulo operation */
    while (len >= NMAX) {
        len -= NMAX;
        n = NMAX / 16;          /* NMAX is divisible by 16 */
        do {
            DO16(buf);          /* 16 sums unrolled */
            buf += 16;
        } while (--n);
        reduce_full(adler);
        reduce_full(sum2);
    }

    /* do remaining bytes (less than NMAX, still just one modulo) */
    if (len) {                  /* avoid modulos if none remaining */
        while (len >= 16) {
            len -= 16;
            DO16(buf);
            buf += 16;
        }
        while (len--) {
            adler += *buf++;
            sum2 += adler;
        }
        reduce_full(adler);
        reduce_full(sum2);
    }

    /* return recombined sums */
    return adler | (sum2 << 16);
}
#  endif
#endif

/* ========================================================================= */
uLong ZEXPORT adler32(adler, buf, len)
    uLong adler;
    const Bytef *buf;
    uInt len;
{
    /* in case user likes doing a byte at a time, keep it fast */
    if (len == 1)
        return adler32_1(adler, buf, len); /* should create a fast tailcall */

    /* initial Adler-32 value (deferred check for len == 1 speed) */
    if (buf == Z_NULL)
        return 1L;

    /* in case short lengths are provided, keep it somewhat fast */
    if (len < MIN_WORK)
        return adler32_common(adler, buf, len);

    return adler32_vec(adler, buf, len);
}

/* ========================================================================= */
local uLong adler32_combine_(adler1, adler2, len2)
    uLong adler1;
    uLong adler2;
    z_off64_t len2;
{
    unsigned long sum1;
    unsigned long sum2;
    unsigned rem;

    /* the derivation of this formula is left as an exercise for the reader */
    rem = (unsigned)(len2 % BASE);
    sum1 = adler1 & 0xffff;
    sum2 = rem * sum1;
    reduce_full(sum2);
    sum1 += (adler2 & 0xffff) + BASE - 1;
    sum2 += ((adler1 >> 16) & 0xffff) + ((adler2 >> 16) & 0xffff) + BASE - rem;
    if (sum1 >= BASE) sum1 -= BASE;
    if (sum1 >= BASE) sum1 -= BASE;
    if (sum2 >= (BASE << 1)) sum2 -= (BASE << 1);
    if (sum2 >= BASE) sum2 -= BASE;
    return sum1 | (sum2 << 16);
}

/* ========================================================================= */
uLong ZEXPORT adler32_combine(adler1, adler2, len2)
    uLong adler1;
    uLong adler2;
    z_off_t len2;
{
    return adler32_combine_(adler1, adler2, len2);
}

uLong ZEXPORT adler32_combine64(adler1, adler2, len2)
    uLong adler1;
    uLong adler2;
    z_off64_t len2;
{
    return adler32_combine_(adler1, adler2, len2);
}
