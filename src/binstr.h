/*
Copyright (C) 2022 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus-3D.

Platanus-3D is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus-3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus-3D; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef __BINSTR_H__
#define __BINSTR_H__
#include <iostream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// 64bit

typedef unsigned long long  u64_t;
const u64_t U64MINUS1 = 0xFFFFFFFFFFFFFFFFull;


////////////////////////////////////////////////////////////////////////////////
// binary string (1char=2bit)

struct binstr_t {
  u64_t               *value;
  const unsigned long long len;  // string length

  //////////////////////////////////////////////////////////////////////////////
  // constructors and destructor

  binstr_t():value(NULL),len(0) {}
  binstr_t(unsigned long long l):len(l) { value=new u64_t[(len+31)/32]; clear(); }
  binstr_t(const binstr_t& rhs):value(NULL),len(rhs.len)
    { value=new u64_t[(len+31)/32];
      u64_t *v=value, *vr=rhs.value;
      for (unsigned long long i=0,n=(len+31)/32; i<n; i++)
        *v++ = *vr++;
    }
  virtual ~binstr_t()
  {
    if(value) {
      delete[] value;
    }
    value = NULL;
  }


  //////////////////////////////////////////////////////////////////////////////
  // utilities
  inline size_t writeTemporaryFile(FILE *fp)
  {
    return fwrite(value, sizeof(u64_t), (len+31)/32, fp);
  }

  inline size_t readTemporaryFile(FILE *fp)
  {
    return fread(value, sizeof(u64_t), (len+31)/32, fp);
  }

  inline void set(const unsigned pos, const unsigned char val)
  {
    value[pos/32] = (value[pos/32] & ~(0x3ull << ((pos%32)*2))) | (static_cast<u64_t>(val) << ((pos%32)*2));
  }

  inline unsigned long long get(unsigned pos)
  {
    return (value[pos/32] >> ((pos%32)*2)) & 3;
  }

  inline void convertToVector(std::vector<u64_t> &vec)
  {
    unsigned long long n = (len + 31) / 32;
    vec.resize(n);
    for (unsigned long long i = 0; i < n; ++i)
      vec[i] = value[i];
  }

  inline void convertFromVector(const std::vector<u64_t> &vec)
  {
    unsigned long long n = (len + 31) / 32;
    for (unsigned long long i = 0; i < n; ++i)
      value[i] = vec[i];
  }

  inline void resize(unsigned long long _len) {
    u64_t *tmp = (_len>0 ? new u64_t[(_len+31)/32] : NULL);
    for (unsigned long long i = 0; i < (_len + 31) / 32; ++i)
        tmp[i] = 0;
    if (value) {
      unsigned long long minLen = this->len < _len ? this->len : _len;
      for (unsigned long long i = 0; i < (minLen + 31) / 32; ++i)
        tmp[i] = value[i];
      delete[] value;
    }
    value = tmp;
    *(const_cast<unsigned long long*>(&(this->len))) = _len;
  }

  inline u64_t to_u64() const {
    u64_t *v=value;
    for (unsigned long long i=1,n=(len+31)/32; i<n; i++) {
      if (*(++v)) {
        return U64MINUS1;
      }
    }
    return *value;
  }

  inline binstr_t& flip() {
    u64_t *v=value;
    for (unsigned long long i=0,n=(len+31)/32; i<n; i++,v++) {
      *v ^= U64MINUS1;
    }
    if (len%32>0) {
      (*(--v)) &= ((0x1ull<<(2*(len%32)))-0x1ull);  // mask
    }
    return *this;
  }

  inline void clear() {
    u64_t *v=value;
    for (unsigned long long i=0,n=(len+31)/32; i<n; i++) {
      *v++ = 0;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // operators
  template<typename T=unsigned>
  inline T operator%(const T &num) const {
    if (!value)
        throw;
    return value[0] % num;
  }

  template<typename T>
  inline binstr_t& operator=(const T rhs) {
    if (!value) {
      resize(1);
    }
    u64_t *v=value;
    for (unsigned long long i=1,n=(len+31)/32; i<n; i++) {
      *(++v) = 0;
    }
    *value = (len>31 ? rhs : (rhs & ((0x1ull<<(2*(len%32)))-0x1ull)));
    return *this;
  }

  inline binstr_t& operator=(const binstr_t& rhs) {
    if (len!=rhs.len || !value) {
      resize(rhs.len);
    }
    u64_t *v=value, *vr=rhs.value;
    for (unsigned long long i=0,n=(len+31)/32; i<n; i++) {
      *v++ = *vr++;
    }
    return *this;
  }
  inline unsigned operator&(unsigned long long num) const {
    return value[0] & num;
  }



#define DEFINE_OPERATOR(B)                                              \
    if (static_cast<long>(len)*2<=static_cast<long>(n) || len==0) {                                           \
      (B).clear();                                                      \
    } else {                                                            \
      unsigned long long i, nb=(len+31)/32, nd=n/64, nm=n%64;                          \
      u64_t *v=(B).value;                                               \
      if (nm==0) {                                                      \
        for (i=0; i<nb-nd; i++) {                                       \
          *v = *(v+nd);                                                 \
          v++;                                                          \
        }                                                               \
      } else {                                                          \
        for (i=0; i<nb-nd-1; i++) {                                     \
          *v = ((*(v+nd))>>nm) | ((*(v+nd+1))<<(64-nm));                \
          v++;                                                          \
        }                                                               \
        *v = (*(v+nd))>>nm;                                             \
        v++;                                                            \
        i++;                                                            \
      }                                                                 \
      for (; i<nb; i++) {                                               \
        *v = 0;                                                         \
        v++;                                                            \
      }                                                                 \
    }

  template <typename TYPE=long>
  inline binstr_t operator>>(TYPE n) const {
    binstr_t b(*this); DEFINE_OPERATOR(b); return b;
  }
  inline binstr_t& operator>>=(long n) {
    DEFINE_OPERATOR(*this); return *this;
  }
#undef DEFINE_OPERATOR

#define DEFINE_OPERATOR(B)                                              \
    if (static_cast<long>(len)*2<=n || len==0) {                                           \
      (B).clear();                                                      \
    } else {                                                            \
      unsigned long long i, nb=(len+31)/32, nd=n/64, nm=n%64;                          \
      u64_t *v=(B).value + nb - 1;                                      \
      if (nm==0) {                                                      \
        for (i=0; i<nb-nd; i++) {                                       \
          *v = *(v-nd);                                                 \
          v--;                                                          \
        }                                                               \
      } else {                                                          \
        for (i=0; i<nb-nd-1; i++) {                                     \
          *v = ((*(v-nd))<<nm) | ((*(v-nd-1))>>(64-nm));                \
          v--;                                                          \
        }                                                               \
        *v = (*(v-nd))<<nm;                                             \
        v--;                                                            \
        i++;                                                            \
      }                                                                 \
      for (; i<nb; i++) {                                               \
        *v = 0;                                                         \
        v--;                                                            \
      }                                                                 \
      if (len%32>0) {                                                   \
        *((B).value+nb-1) &= ((0x1ull<<(2*(len%32)))-0x1ull);  /* mask */ \
      }                                                                 \
    }

  template <typename TYPE=long>
  inline binstr_t operator<<(TYPE n) const {
    binstr_t b(*this); DEFINE_OPERATOR(b); return b;
  }
  inline binstr_t& operator<<=(long n) {
    DEFINE_OPERATOR(*this); return *this;
  }
#undef DEFINE_OPERATOR

  inline bool operator==(const binstr_t &rhs) const {
    if (len != rhs.len) return false;
    u64_t *v=value, *vr=rhs.value;
    for (u64_t i = 0, n = (len+31)/32; i < n; ++i) {
      if (*v++ != *vr++) {
        return false;
      }
    }
    return true;
  }

  inline bool operator<(const binstr_t& rhs) const {
    unsigned long long i=(len+31)/32-1;
    u64_t *l=value, *r=rhs.value;
    while (i && *(l+i)==*(r+i)) { i--; }
    return (*(l+i)<*(r+i));
  };
  // hash function for unordered_map
  struct hasher
  {
    inline size_t operator()(binstr_t b) const { return (size_t)*(b.value); }
  };
  static void revCom_u64(unsigned long long &t) {
    t ^= UINT64_MAX;
    t = ((t&0x3333333333333333ull)<<2) | ((t&0xCCCCCCCCCCCCCCCCull)>>2);
    t = ((t&0x0F0F0F0F0F0F0F0Full)<<4) | ((t&0xF0F0F0F0F0F0F0F0ull)>>4);
    t = ((t&0x00FF00FF00FF00FFull)<<8) | ((t&0xFF00FF00FF00FF00ull)>>8);
    t = ((t&0x0000FFFF0000FFFFull)<<16) | ((t&0xFFFF0000FFFF0000ull)>>16);
    t = (t<<32) | (t>>32);
  }

};







////////////////////////////////////////////////////////////////////////////////
// binary string (1char=2bit)

struct BinstrBase {
  u64_t               *value;
  const unsigned long long len;  // string length

  //////////////////////////////////////////////////////////////////////////////
  // constructors and destructor

  BinstrBase():value(NULL),len(0) {}
  BinstrBase(unsigned long long l):value(NULL), len(l) {}
  BinstrBase(const BinstrBase& rhs):value(NULL),len(rhs.len)
  {
  }
  virtual ~BinstrBase()
  {
  }


  //////////////////////////////////////////////////////////////////////////////
  // utilities
  inline size_t writeTemporaryFile(FILE *fp)
  {
    return fwrite(value, sizeof(u64_t), (len+31)/32, fp);
  }

  inline size_t readTemporaryFile(FILE *fp)
  {
    return fread(value, sizeof(u64_t), (len+31)/32, fp);
  }

  inline void set(const unsigned pos, const unsigned char val)
  {
    value[pos/32] = (value[pos/32] & ~(0x3ull << ((pos%32)*2))) | (static_cast<u64_t>(val) << ((pos%32)*2));
  }

  inline unsigned long long get(unsigned pos)
  {
    return (value[pos/32] >> ((pos%32)*2)) & 3;
  }

  inline void convertToVector(std::vector<u64_t> &vec)
  {
    unsigned long long n = (len + 31) / 32;
    vec.resize(n);
    for (unsigned long long i = 0; i < n; ++i)
      vec[i] = value[i];
  }

  inline void convertFromVector(const std::vector<u64_t> &vec)
  {
    unsigned long long n = (len + 31) / 32;
    for (unsigned long long i = 0; i < n; ++i)
      value[i] = vec[i];
  }

  inline void clear() {
    u64_t *v=value;
    for (unsigned long long i=0,n=(len+31)/32; i<n; i++) {
      *v++ = 0;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // operators
  template<typename T=unsigned>
  inline T operator%(const T &num) const {
    if (!value)
        throw;
    return value[0] % num;
  }

  inline virtual BinstrBase& operator=(const BinstrBase& rhs) {
    *(const_cast<unsigned long long*>(&(this->len))) = rhs.len;
    u64_t *v=value;
    const u64_t *vr=rhs.value;
    for (unsigned long long i=0,n=(len+31)/32; i<n; i++) {
      *v++ = *vr++;
    }
    return *this;
  }
  inline unsigned operator&(unsigned long long num) const {
    return *value & num;
  }

#define DEFINE_OPERATOR(B)                                              \
    if (static_cast<long>(len)*2<=static_cast<long>(n) || len==0) {                                           \
      (B).clear();                                                      \
    } else {                                                            \
      unsigned long long i, nb=(len+31)/32, nd=n/64, nm=n%64;                          \
      u64_t *v=(B).value;                                               \
      if (nm==0) {                                                      \
        for (i=0; i<nb-nd; i++) {                                       \
          *v = *(v+nd);                                                 \
          v++;                                                          \
        }                                                               \
      } else {                                                          \
        for (i=0; i<nb-nd-1; i++) {                                     \
          *v = ((*(v+nd))>>nm) | ((*(v+nd+1))<<(64-nm));                \
          v++;                                                          \
        }                                                               \
        *v = (*(v+nd))>>nm;                                             \
        v++;                                                            \
        i++;                                                            \
      }                                                                 \
      for (; i<nb; i++) {                                               \
        *v = 0;                                                         \
        v++;                                                            \
      }                                                                 \
    }

  template <typename TYPE=long>
  inline BinstrBase operator>>(TYPE n) const {
    BinstrBase b(*this);DEFINE_OPERATOR(b); return b;
  }
  inline virtual BinstrBase& operator>>=(long n) {
    DEFINE_OPERATOR(*this); return *this;
  }
#undef DEFINE_OPERATOR

#define DEFINE_OPERATOR(B)                                              \
    if (static_cast<long>(len)*2<=n || len==0) {                                           \
      (B).clear();                                                      \
    } else {                                                            \
      unsigned long long i, nb=(len+31)/32, nd=n/64, nm=n%64;                          \
      u64_t *v=(B).value + nb - 1;                                      \
      if (nm==0) {                                                      \
        for (i=0; i<nb-nd; i++) {                                       \
          *v = *(v-nd);                                                 \
          v--;                                                          \
        }                                                               \
      } else {                                                          \
        for (i=0; i<nb-nd-1; i++) {                                     \
          *v = ((*(v-nd))<<nm) | ((*(v-nd-1))>>(64-nm));                \
          v--;                                                          \
        }                                                               \
        *v = (*(v-nd))<<nm;                                             \
        v--;                                                            \
        i++;                                                            \
      }                                                                 \
      for (; i<nb; i++) {                                               \
        *v = 0;                                                         \
        v--;                                                            \
      }                                                                 \
      if (len%32>0) {                                                   \
        *((B).value+nb-1) &= ((0x1ull<<(2*(len%32)))-0x1ull);  /* mask */ \
      }                                                                 \
    }

  template <typename TYPE=long>
  inline BinstrBase operator<<(TYPE n) const {
    BinstrBase b(*this); DEFINE_OPERATOR(b); return b;
  }
  inline virtual BinstrBase& operator<<=(long n) {
    DEFINE_OPERATOR(*this); return *this;
  }
#undef DEFINE_OPERATOR

  inline bool operator==(const BinstrBase &rhs) const {
    if (len != rhs.len) return false;
    const u64_t *v=value;
    const u64_t *vr=rhs.value;
    for (u64_t i = 0, n = (len+31)/32; i < n; ++i) {
      if (*v++ != *vr++) {
        return false;
      }
    }
    return true;
  }

  inline bool operator<(const BinstrBase& rhs) const {
    unsigned long long i=(len+31)/32-1;
    const u64_t *l=value;
    const u64_t *r=rhs.value;
    while (i && *(l+i)==*(r+i)) { i--; }
    return (*(l+i)<*(r+i));
  };

  static void revCom_u64(unsigned long long &t) {
    t ^= UINT64_MAX;
    t = ((t&0x3333333333333333ull)<<2) | ((t&0xCCCCCCCCCCCCCCCCull)>>2);
    t = ((t&0x0F0F0F0F0F0F0F0Full)<<4) | ((t&0xF0F0F0F0F0F0F0F0ull)>>4);
    t = ((t&0x00FF00FF00FF00FFull)<<8) | ((t&0xFF00FF00FF00FF00ull)>>8);
    t = ((t&0x0000FFFF0000FFFFull)<<16) | ((t&0xFFFF0000FFFF0000ull)>>16);
    t = (t<<32) | (t>>32);
  }
  struct hasher
  {
    inline size_t operator()(BinstrBase b) const { return (size_t)*(b.value); }
  };
};



struct Binstr63 : public BinstrBase
{
  u64_t entity[2];

  Binstr63():BinstrBase(),entity() {value = entity; clear();}
  Binstr63(const unsigned long long l):BinstrBase(l),entity() {value = entity; clear();}
  Binstr63(const Binstr63& rhs):BinstrBase(rhs),entity()
  {
    value = entity;
    u64_t *v=value;
    const u64_t *vr=rhs.value;
    for (unsigned long long i=0; i<2; i++) {
      *v++ = *vr++;
    }
  }


  virtual ~Binstr63()
  {
  }


  inline Binstr63 &operator=(const Binstr63 &rhs)
  {
    BinstrBase::operator=(rhs);
    return *this;
  }

  template <typename TYPE=long>
  inline Binstr63 operator>>(TYPE n) const {
    Binstr63 b(*this); b >>= n; return b;
  }

  inline Binstr63 &operator>>=(long n)
  {
    BinstrBase::operator>>=(n);
    return *this;
  }

  template <typename TYPE=long>
  inline Binstr63 operator<<(TYPE n) const {
    Binstr63 b(*this); b <<= n; return b;
  }

  inline Binstr63 &operator<<=(long n)
  {
    BinstrBase::operator<<=(n);
    return *this;
  }


  struct hasher
  {
    inline size_t operator()(Binstr63 b) const { return (size_t)(b.entity[0]); }
  };

};



struct Binstr95 : public BinstrBase
{
  u64_t entity[3];

  Binstr95():BinstrBase(),entity() {value = entity; clear();}
  Binstr95(const unsigned long long l):BinstrBase(l),entity() {value = entity; clear();}
  Binstr95(const Binstr95& rhs):BinstrBase(rhs),entity()
  {
    value = entity;
    u64_t *v=value;
    const u64_t *vr=rhs.value;
    for (unsigned long long i=0; i<3; i++) {
      *v++ = *vr++;
    }
  }


  virtual ~Binstr95()
  {
  }


  inline Binstr95 &operator=(const Binstr95 &rhs)
  {
    BinstrBase::operator=(rhs);
    return *this;
  }

  template <typename TYPE=long>
  inline Binstr95 operator>>(TYPE n) const {
    Binstr95 b(*this); b >>= n; return b;
  }

  inline Binstr95 &operator>>=(long n)
  {
    BinstrBase::operator>>=(n);
    return *this;
  }

  template <typename TYPE=long>
  inline Binstr95 operator<<(TYPE n) const {
    Binstr95 b(*this); b <<= n; return b;
  }

  inline Binstr95 &operator<<=(long n)
  {
    BinstrBase::operator<<=(n);
    return *this;
  }


  struct hasher
  {
    inline size_t operator()(Binstr95 b) const { return (size_t)(b.entity[0]); }
  };

};



struct Binstr127 : public BinstrBase
{
  u64_t entity[4];

  Binstr127():BinstrBase(),entity() {value = entity; clear();}
  Binstr127(const unsigned long long l):BinstrBase(l),entity() {value = entity; clear();}
  Binstr127(const Binstr127& rhs):BinstrBase(rhs),entity()
  {
    value = entity;
    u64_t *v=value;
    const u64_t *vr=rhs.value;
    for (unsigned long long i=0; i<4; i++) {
      *v++ = *vr++;
    }
  }


  virtual ~Binstr127()
  {
  }


  inline Binstr127 &operator=(const Binstr127 &rhs)
  {
    BinstrBase::operator=(rhs);
    return *this;
  }

  template <typename TYPE=long>
  inline Binstr127 operator>>(TYPE n) const {
    Binstr127 b(*this); b >>= n; return b;
  }

  inline Binstr127 &operator>>=(long n)
  {
    BinstrBase::operator>>=(n);
    return *this;
  }

  template <typename TYPE=long>
  inline Binstr127 operator<<(TYPE n) const {
    Binstr127 b(*this); b <<= n; return b;
  }

  inline Binstr127 &operator<<=(long n)
  {
    BinstrBase::operator<<=(n);
    return *this;
  }


  struct hasher
  {
    inline size_t operator()(Binstr127 b) const { return (size_t)(b.entity[0]); }
  };

};




struct Binstr159 : public BinstrBase
{
  u64_t entity[5];

  Binstr159():BinstrBase(),entity() {value = entity; clear();}
  Binstr159(const unsigned long long l):BinstrBase(l),entity() {value = entity; clear();}
  Binstr159(const Binstr159& rhs):BinstrBase(rhs),entity()
  {
    value = entity;
    u64_t *v=value;
    const u64_t *vr=rhs.value;
    for (unsigned long long i=0; i<5; i++) {
      *v++ = *vr++;
    }
  }


  virtual ~Binstr159()
  {
  }


  inline Binstr159 &operator=(const Binstr159 &rhs)
  {
    BinstrBase::operator=(rhs);
    return *this;
  }

  template <typename TYPE=long>
  inline Binstr159 operator>>(TYPE n) const {
    Binstr159 b(*this); b >>= n; return b;
  }

  inline Binstr159 &operator>>=(long n)
  {
    BinstrBase::operator>>=(n);
    return *this;
  }

  template <typename TYPE=long>
  inline Binstr159 operator<<(TYPE n) const {
    Binstr159 b(*this); b <<= n; return b;
  }

  inline Binstr159 &operator<<=(long n)
  {
    BinstrBase::operator<<=(n);
    return *this;
  }


  struct hasher
  {
    inline size_t operator()(Binstr159 b) const { return (size_t)(b.entity[0]); }
  };

};

#endif
