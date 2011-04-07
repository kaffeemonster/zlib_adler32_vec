/*
 * adler32.c -- compute the Adler-32 checksum of a data stream
 *   arm implementation
 * Copyright (C) 1995-2007 Mark Adler
 * Copyright (C) 2009-2011 Jan Seiffert
 * For conditions of distribution and use, see copyright notice in zlib.h
 */

/* @(#) $Id$ */

#if defined(__ARM_NEON__)
// TODO: need byte order define
/*
 * Big endian NEON qwords are kind of broken.
 * They are big endian within the dwords, but WRONG
 * (really??) way round between lo and hi.
 * Creating some kind of PDP11 middle endian.
 *
 * This is madness and unsupportable. For this reason
 * GCC wants to disable qword endian specific patterns.
 * We would need a Preprocessor define which endian we
 * have to disable this code.
 */
#  include <arm_neon.h>

#  define SOVUCQ sizeof(uint8x16_t)
#  define SOVUC sizeof(uint8x8_t)
/* since we do not have the 64bit psadbw sum, we could prop. do a little more */
#  define VNMAX (6*NMAX)
#  define HAVE_ADLER32_VEC
#  define MIN_WORK 32

/* ========================================================================= */
local inline uint8x16_t neon_simple_alignq(uint8x16_t a, uint8x16_t b, unsigned amount)
{
    switch(amount % SOVUCQ)
    {
    case  0: return a;
    case  1: return vextq_u8(a, b,  1);
    case  2: return vextq_u8(a, b,  2);
    case  3: return vextq_u8(a, b,  3);
    case  4: return vextq_u8(a, b,  4);
    case  5: return vextq_u8(a, b,  5);
    case  6: return vextq_u8(a, b,  6);
    case  7: return vextq_u8(a, b,  7);
    case  8: return vextq_u8(a, b,  8);
    case  9: return vextq_u8(a, b,  9);
    case 10: return vextq_u8(a, b, 10);
    case 11: return vextq_u8(a, b, 11);
    case 12: return vextq_u8(a, b, 12);
    case 13: return vextq_u8(a, b, 13);
    case 14: return vextq_u8(a, b, 14);
    case 15: return vextq_u8(a, b, 15);
    }
    return b;
}

/* ========================================================================= */
local inline uint32x4_t vector_reduce(uint32x4_t x)
{
    uint32x4_t y;

    y = vshlq_n_u32(x, 16);
    x = vshrq_n_u32(x, 16);
    y = vshrq_n_u32(y, 16);
    y = vsubq_u32(y, x);
    x = vaddq_u32(y, vshlq_n_u32(x, 4));
    return x;
}

/* ========================================================================= */
local noinline uLong adler32_vec(adler, buf, len)
    uLong adler;
    const Bytef *buf;
    uInt len;
{
    uint32x4_t v0_32 = (uint32x4_t){0,0,0,0};
    uint8x16_t    v0 = (uint8x16_t)v0_32;
    uint8x16_t vord, vord_a;
    uint32x4_t vs1, vs2;
    uint32x2_t v_tsum;
    uint8x16_t in16;
    uint32_t s1, s2;
    unsigned k;

    s1 = adler & 0xffff;
    s2 = (adler >> 16) & 0xffff;

// TODO: big endian mask is prop. wrong
    if (host_is_bigendian())
        vord = (uint8x16_t){16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
    else
        vord = (uint8x16_t){1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

    if (likely(len >= 2*SOVUCQ)) {
        unsigned f, n;

        /*
         * Add stuff to achieve alignment
         */
        /* align hard down */
        f = (unsigned) ALIGN_DOWN_DIFF(buf, SOVUCQ);
        n = SOVUCQ - f;
        buf = (const unsigned char *)ALIGN_DOWN(buf, SOVUCQ);

        /* add n times s1 to s2 for start round */
        s2 += s1 * n;

        /* set sums 0 */
        vs1 = v0_32;
        vs2 = v0_32;
        /*
         * the accumulation of s1 for every round grows very fast
         * (quadratic?), even if we accumulate in 4 dwords, more
         * rounds means nonlinear growth.
         * We already split it out of s2, normaly it would be in
         * s2 times 16... and even grow faster.
         * Thanks to this split and vector reduction, we can stay
         * longer in the loops. But we have to prepare for the worst
         * (all 0xff), only do 6 times the work.
         * (we could prop. stay a little longer since we have 4 sums,
         * not 2 like on x86).
         */
        k = len < VNMAX ? (unsigned)len : VNMAX;
        len -= k;
        /* insert scalar start somewhere */
        vs1 = vsetq_lane_u32(s1, vs1, 0);
        vs2 = vsetq_lane_u32(s2, vs2, 0);

        /* get input data */
        in16 = *(const uint8x16_t *)buf;
        /* mask out excess data */
        if(host_is_bigendian()) {
            in16 = neon_simple_alignq(v0, in16, n);
            vord_a = neon_simple_alignq(v0, vord, n);
        } else {
            in16 = neon_simple_alignq(in16, v0, f);
            vord_a = neon_simple_alignq(vord, v0, f);
        }

        /* pairwise add bytes and long, pairwise add word long acc */
        vs1 = vpadalq_u16(vs1, vpaddlq_u8(in16));
        /* apply order, add words, pairwise add word long acc */
        vs2 = vpadalq_u16(vs2,
                vmlal_u8(
                    vmull_u8(vget_low_u8(in16), vget_low_u8(vord_a)),
                    vget_high_u8(in16), vget_high_u8(vord_a)
                    )
                );

        buf += SOVUCQ;
        k -= n;

        if (likely(k >= SOVUCQ)) do {
            uint32x4_t vs1_r = v0_32;
            do {
                /* add vs1 for this round */
                vs1_r = vaddq_u32(vs1_r, vs1);

                /* get input data */
                in16 = *(const uint8x16_t *)buf;

// TODO: make work in inner loop more tight
                /*
                 * decompose partial sums, so we do less instructions and
                 * build loops around it to do acc and so on only from time
                 * to time.
                 * This is hard with NEON, because the instruction are nice:
                 * we have the stuff in widening and with acc (practicaly
                 * for free...)
                 */
                /* pairwise add bytes and long, pairwise add word long acc */
                vs1 = vpadalq_u16(vs1, vpaddlq_u8(in16));
                /* apply order, add words, pairwise add word long acc */
                vs2 = vpadalq_u16(vs2,
                        vmlal_u8(
                            vmull_u8(vget_low_u8(in16), vget_low_u8(vord)),
                            vget_high_u8(in16), vget_high_u8(vord)
                            )
                        );

                buf += SOVUCQ;
                k -= SOVUCQ;
            } while (k >= SOVUCQ);
            /* reduce vs1 round sum before multiplying by 16 */
            vs1_r = vector_reduce(vs1_r);
            /* add vs1 for this round (16 times) */
            /* they have shift right and accummulate, where is shift left and acc?? */
            vs2 = vaddq_u32(vs2, vshlq_n_u32(vs1_r, 4));
            /* reduce both vectors to something within 16 bit */
            vs2 = vector_reduce(vs2);
            vs1 = vector_reduce(vs1);
            len += k;
            k = len < VNMAX ? (unsigned) len : VNMAX;
            len -= k;
        } while (likely(k >= SOVUC));

        if (likely(k)) {
            /*
             * handle trailer
             */
            f = SOVUCQ - k;
            /* add k times vs1 for this trailer */
            vs2 = vmlaq_u32(vs2, vs1, vdupq_n_u32(k));

            /* get input data */
            in16 = *(const uint8x16_t *)buf;
            /* masks out bad data */
            if(host_is_bigendian())
                in16 = neon_simple_alignq(in16, v0, f);
            else
                in16 = neon_simple_alignq(v0, in16, k);

            /* pairwise add bytes and long, pairwise add word long acc */
            vs1 = vpadalq_u16(vs1, vpaddlq_u8(in16));
            /* apply order, add words, pairwise add word long acc */
            vs2 = vpadalq_u16(vs2,
                    vmlal_u8(
                        vmull_u8(vget_low_u8(in16), vget_low_u8(vord)),
                        vget_high_u8(in16), vget_high_u8(vord)
                        )
                    );

            buf += k;
            k -= k;
        }

        /* add horizontal */
        v_tsum = vpadd_u32(vget_high_u32(vs1), vget_low_u32(vs1));
        v_tsum = vpadd_u32(v_tsum, v_tsum);
        s1 = vget_lane_u32(v_tsum, 0);
        v_tsum = vpadd_u32(vget_high_u32(vs2), vget_low_u32(vs2));
        v_tsum = vpadd_u32(v_tsum, v_tsum);
        s2 = vget_lane_u32(v_tsum, 0);
    }

    if (unlikely(len)) do {
        s1 += *buf++;
        s2 += s1;
    } while (--len);
    reduce_x(s1);
    reduce_x(s2);

    return (s2 << 16) | s1;
}

/* inline asm, so only on GCC (or compatible) && ARM v6 or better */
#elif defined(__GNUC__) && ( \
        defined(__ARM_ARCH_6__)  || defined(__ARM_ARCH_6J__)  || \
        defined(__ARM_ARCH_6Z__) || defined(__ARM_ARCH_6ZK__) || \
        defined(__ARM_ARCH_7A__) \
      )
#  define SOU32 (sizeof(unsigned int))
#  define HAVE_ADLER32_VEC
#  define MIN_WORK 16

/* ========================================================================= */
local noinline uLong adler32_vec(adler, buf, len)
    uLong adler;
    const Bytef *buf;
    uInt len;
{
    unsigned int s1, s2;
    unsigned int k;

    s1 = adler & 0xffff;
    s2 = (adler >> 16) & 0xffff;

    k    = ALIGN_DIFF(buf, SOU32);
    len -= k;
    if (k) do {
        s1 += *buf++;
        s2 += s1;
    } while (--k);

    if (likely(len >= 4 * SOU32)) {
        unsigned int vs1 = s1, vs2 = s2;
        unsigned int order_lo, order_hi;

// TODO: byte order?
        if(host_is_bigendian()) {
            order_lo = 0x00030001;
            order_hi = 0x00040002;
        } else {
            order_lo = 0x00020004;
            order_hi = 0x00010003;
        }
// TODO: we could go over NMAX, since we have split the vs2 sum
        /* something around (NMAX+(NMAX/3)+302) */
        k = len < NMAX ? len : NMAX;
        len -= k;

        do {
            unsigned int vs1_r = 0;
            do {
                unsigned int t21, t22, in;

                /* get input data */
                in = *(const unsigned int *)buf;

                /* add vs1 for this round */
                vs1_r += vs1;

                /* add horizontal and acc */
                asm ("usada8 %0, %1, %2, %3" : "=r" (vs1) : "r" (in), "r" (0), "r" (vs1));
                /* widen bytes to words, apply order, add and acc */
                asm ("uxtb16 %0, %1" : "=r" (t21) : "r" (in));
                asm ("uxtb16 %0, %1, ror #8" : "=r" (t22) : "r" (in));
// TODO: instruction result latency
                /*
                 * The same problem like the classic serial sum:
                 * Chip maker sell us 1-cycle instructions, but that is not the
                 * whole story. Nearly all 1-cycle chips are pipelined, so
                 * you can get one result per cycle, but only if _they_ (plural)
                 * are independent.
                 * If you are depending on the result of an preciding instruction,
                 * in the worst case you hit the instruction latency which is worst
                 * case >= pipeline length. On the other hand there are result-fast-paths.
                 * This could all be a wash with the classic sum (4 * 2 instructions,
                 * + dependence), since smald is:
                 * - 2 cycle issue
                 * - needs the acc in pipeline step E1, instead of E2
                 * But the Cortex has a fastpath for acc.
                 * I don't know.
                 * We can not even unroll, we would need 4 order vars, return ENOREGISTER.
                 */
                asm ("smlad %0, %1, %2, %3" : "=r" (vs2) : "r" (t21) , "r" (order_lo), "r" (vs2));
                asm ("smlad %0, %1, %2, %3" : "=r" (vs2) : "r" (t22) , "r" (order_hi), "r" (vs2));

                buf += SOU32;
                k -= SOU32;
            } while (k >= SOU32);
            /* reduce vs1 round sum before multiplying by 4 */
            reduce(vs1_r);
            /* add vs1 for this round (4 times) */
            vs2 += vs1_r * 4;
            /* reduce both sums to something within 16 bit */
            reduce(vs2);
            reduce(vs1);
            len += k;
            k = len < NMAX ? len : NMAX;
            len -= k;
        } while (likely(k >= 4 * SOU32));
        len += k;
        s1 = vs1;
        s2 = vs2;
    }

    if (unlikely(len)) do {
        s1 += *buf++;
        s2 += s1;
    } while (--len);
    /* at this point we should no have so big s1 & s2 */
    reduce_x(s1);
    reduce_x(s2);

    return (s2 << 16) | s1;
}
#endif
