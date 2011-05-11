/* slhash.c -- slide the hash table during fill_window()
 * Copyright (C) 1995-2010 Jean-loup Gailly and Mark Adler
 * Copyright (C) 2011 Jan Seiffert
 * For conditions of distribution and use, see copyright notice in zlib.h
 */

/* NOTE:
 * We do not precheck the length or wsize for small values because
 * we assume a minimum len of 256 (for MEM_LEVEL 1) and a minimum wsize
 * of 256 for windowBits 8
 */

#if defined(__ARM_NEON__) && defined(__ARMEL__)
/*
 * Big endian NEON qwords are kind of broken.
 * They are big endian within the dwords, but WRONG
 * (really??) way round between lo and hi.
 * Creating some kind of PDP11 middle endian.
 *
 * This is madness and unsupportable. For this reason
 * GCC wants to disable qword endian specific patterns.
 */
#  include <arm_neon.h>

#  define SOVUSQ sizeof(uint16x8_t)
#  define SOVUS sizeof(uint16x4_t)
#  define HAVE_SLHASH_VEC

local void update_hoffset_l(Posf *p, uInt wsize, unsigned n)
{
    unsigned int i, j;
    uint32x4_t vwsize;

    vwsize = vdupq_n_u32(wsize);

    if (ALIGN_DOWN_DIFF(p, SOVUSQ)) {
        uint32x4_t in = vmovl_u16(*(uint16x4_t *)p);
        in  = vqsubq_u32(in, vwsize);
        *(uint16x4_t *)p = vmovn_u32(in);
        p += SOVUS/sizeof(*p);
        n -= SOVUS/sizeof(*p);
    }
    i  = n / (SOVUSQ/sizeof(*p));
    n %= SOVUSQ/sizeof(*p);

    do {
        uint16x8_t in8 = *(uint16x8_t *)p;
        uint32x4_t inl = vmovl_u16(vget_low_u16(in8));
        uint32x4_t inh = vmovl_u16(vget_high_u16(in8));
        inl = vqsubq_u32(inl, vwsize);
        inh = vqsubq_u32(inh, vwsize);
        in8 = vcombine_u16(vmovn_u32(inl), vmovn_u32(inh));
        *(uint16x8_t *)p = in8;
        p += SOVUSQ/sizeof(*p);
    } while (--i);

    if (n >= SOVUS/sizeof(*p)) {
        uint32x4_t in = vmovl_u16(*(uint16x4_t *)p);
        in  = vqsubq_u32(in, vwsize);
        *(uint16x4_t *)p = vmovn_u32(in);
        p += SOVUS/sizeof(*p);
        n -= SOVUS/sizeof(*p);
    }
    if (unlikely(n)) do {
        register unsigned m = *p;
        *p++ = (Pos)(m >= wsize ? m-wsize : NIL);
    } while (--n);
}

local void update_hoffset(Posf *p, uInt wsize, unsigned n)
{
    unsigned int i;
    uint16x8_t vwsize;

    i  = ALIGN_DIFF(p, SOVUS)/sizeof(Pos);
    n -= i;
    if (unlikely(i)) do {
        register unsigned m = *p;
        *p++ = (Pos)(m >= wsize ? m-wsize : NIL);
    } while (--i);

    if(unlikely(wsize > (1<<16)-1)) {
        update_hoffset_l(p, wsize, n);
        return;
    }
    vwsize = vdupq_n_u16(wsize);

    if (ALIGN_DOWN_DIFF(p, SOVUSQ)) {
        uint16x4_t in4 = *(uint16x4_t *)p;
        in4 = vqsub_u16(in4, vget_low_u16(vwsize));
        *(uint16x4_t *)p = in4;
        p += SOVUS/sizeof(*p);
        n -= SOVUS/sizeof(*p);
    }
    i  = n / (SOVUSQ/sizeof(*p));
    n %= SOVUSQ/sizeof(*p);
    do {
        *(uint16x8_t *)p = vqsubq_u16(*(uint16x8_t *)p, vwsize);
        p += SOVUSQ/sizeof(*p);
    } while (--i);

    if (n >= SOVUS/sizeof(*p)) {
        uint16x4_t in4 = *(uint16x4_t *)p;
        in4 = vqsub_u16(in4, vget_low_u16(vwsize));
        *(uint16x4_t *)p = in4;
        p += SOVUS/sizeof(*p);
        n -= SOVUS/sizeof(*p);
    }
    if (unlikely(n)) do {
        register unsigned m = *p;
        *p++ = (Pos)(m >= wsize ? m-wsize : NIL);
    } while (--n);
}
#endif
