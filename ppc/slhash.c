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

#if defined(__ALTIVEC__) && defined(__GNUC__)
#  define HAVE_SLHASH_VEC
#  include <altivec.h>
#  define SOVUS (sizeof(vector unsigned short))

local inline vector int vector_splat_scalar_s32(int x)
{
    vector int val = vec_lde(0, &x);
    vector unsigned char vperm;

    vperm = (vector unsigned char)vec_splat((vector int)vec_lvsl(0, &x), 0);
    return vec_perm(val, val, vperm);
}

local inline vector unsigned short vector_splat_scalar_u16(unsigned short x)
{
    vector unsigned short val = vec_lde(0, &x);
    vector unsigned char vperm;

    vperm = (vector unsigned char)vec_splat((vector unsigned short)vec_lvsl(0, &x), 0);
    return vec_perm(val, val, vperm);
}

local void update_hoffset_l(Posf *p, uInt wsize, unsigned n)
{
    vector unsigned short vin;
    vector unsigned short v0;
    vector int vwsize, vinl, vinh;
    unsigned int f, i;

    v0 = vec_splat_u16(0);
    vwsize = vector_splat_scalar_s32((int)wsize);
    /* achive alignment */
    f = (unsigned) ALIGN_DOWN_DIFF(p, SOVUS);
    if (f) {
        vector unsigned char vperml, vpermr;
        vector unsigned short vin_o;
            f  = (SOVUS - f)/sizeof(*p);
            n -= f;
        vperml = vec_lvsl(0, p);
        vpermr = vec_lvsr(0, p);
        /* align hard down */
             p = (unsigned short *)ALIGN_DOWN(p, SOVUS);
           vin = vec_ldl(0, p);
         vin_o = vec_perm(vin, vin, vperml);
           vin = vec_perm(vin, v0, vperml);
          vinl = (vector int)vec_mergel(v0, vin);
          vinh = (vector int)vec_mergeh(v0, vin);
          vinl = vec_sub(vinl, vwsize);
          vinh = vec_sub(vinh, vwsize);
           vin = vec_packsu(vinh, vinl);
           vin = vec_perm(vin_o, vin, vpermr);
        /*
         * We write it out again element by element
         * because writing out the whole vector
         * may race against the owner of the other
         * data.
         * The proper fix for this is aligning the
         * arrays at malloc time.
         */
        do {
            vec_ste(vin, SOVUS-(f*sizeof(*p)), p);
        } while (--f);
        p += SOVUS/sizeof(*p);
    }

    i  = n / (SOVUS/sizeof(*p));
    n %= SOVUS/sizeof(*p);

    do {
         vin = vec_ldl(0, p);
        vinl = (vector int)vec_mergel(v0, vin);
        vinh = (vector int)vec_mergeh(v0, vin);
        vinl = vec_sub(vinl, vwsize);
        vinh = vec_sub(vinh, vwsize);
         vin = vec_packsu(vinh, vinl);
        vec_stl(vin, 0, p);
          p += SOVUS/sizeof(*p);
    } while (--i);

    if (n) {
        vector unsigned char vperml, vpermr;
        vector unsigned short vin_o;
           vin = vec_ldl(0, p);
        vperml = vec_lvsl(n*sizeof(*p), p);
        vpermr = vec_lvsr(n*sizeof(*p), p);
         vin_o = vec_perm(vin, vin, vperml);
           vin = vec_perm(v0, vin, vperml);
          vinl = (vector int)vec_mergel(v0, vin);
          vinh = (vector int)vec_mergeh(v0, vin);
          vinl = vec_sub(vinl, vwsize);
          vinh = vec_sub(vinh, vwsize);
           vin = vec_packsu(vinh, vinl);
           vin = vec_perm(vin, vin_o, vpermr);
             i = 0;
        do {
            vec_ste(vin, i*sizeof(*p), p);
            i++;
        } while (--n);
    }
}

local void update_hoffset(Posf *p, uInt wsize, unsigned n)
{
    vector unsigned short vin;
    vector unsigned short v0;
    vector unsigned short vwsize;
    unsigned int f, i;
    unsigned int block_num;

    block_num = DIV_ROUNDUP(n, 512); /* 32 block size * 16 */
    f  = 512;
    f |= block_num >= 256 ? 0 : block_num << 16;
    vec_dst(p, f, 2);

    if (unlikely(wsize > (1<<16)-1)) {
        update_hoffset_l(p, wsize, n);
        return;
    }

        v0 = vec_splat_u16(0);
    vwsize = vector_splat_scalar_u16((unsigned short)wsize);
    /* achive alignment */
    f = (unsigned) ALIGN_DOWN_DIFF(p, SOVUS);
    if (f) {
        vector unsigned char vperml, vpermr;
        vector unsigned short vin_o;
            f  = (SOVUS-f)/sizeof(*p);
            n -= f;
        vperml = vec_lvsl(0, p);
        vpermr = vec_lvsr(0, p);
        /* align hard down */
             p = (unsigned short *)ALIGN_DOWN(p, SOVUS);
           vin = vec_ldl(0, p);
         vin_o = vec_perm(vin, vin, vperml);
           vin = vec_perm(vin, v0, vperml);
           vin = vec_subs(vin, vwsize);
           vin = vec_perm(vin_o, vin, vpermr);
        do {
            vec_ste(vin, SOVUS-(f*sizeof(*p)), p);
        } while (--f);
        p += SOVUS/sizeof(*p);
    }

    i  = n / (SOVUS/sizeof(*p));
    n %= SOVUS/sizeof(*p);

    do {
        vin = vec_ldl(0, p);
        vin = vec_subs(vin, vwsize);
        vec_stl(vin, 0, p);
         p += SOVUS/sizeof(*p);
    } while (--i);

    if (n) {
        vector unsigned char vperml, vpermr;
        vector unsigned short vin_o;
           vin = vec_ldl(0, p);
        vperml = vec_lvsl(n*sizeof(*p), p);
        vpermr = vec_lvsr(n*sizeof(*p), p);
         vin_o = vec_perm(vin, vin, vperml);
           vin = vec_perm(v0, vin, vperml);
           vin = vec_subs(vin, vwsize);
           vin = vec_perm(vin, vin_o, vpermr);
             i = 0;
        do {
            vec_ste(vin, i*sizeof(*p), p);
            i++;
        } while (--n);
    }
}
#endif