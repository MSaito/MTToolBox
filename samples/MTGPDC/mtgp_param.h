#ifndef MTGP_PARAM_H
#define MTGP_PARAM_H
/**
 * @file mtgp_param.h
 *
 * @brief a class for the parameter of mtgps.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (c) 2010 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 * Copyright (c) 2011 Mutsuo Saito, Makoto Matsumoto, Hiroshima
 * University and University of Tokyo. All rights reserved.
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */

#include <stdint.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <MTToolBox/util.hpp>

#if defined (USE_SHA)
#include <openssl/sha.h>
#endif
namespace mtgp {
    template<typename T> class mtgp_param {
    public:
        mtgp_param() {
            for (int i = 0; i < 4; i++) {
                tbl[i] = 0;
                tmp_tbl[i] = 0;
            }
            for (int i = 0; i < 16; i++) {
                p[i] = 0;
                tp[i] = 0;
            }
        }
        mtgp_param(const mtgp_param<T>& that) {
            mexp = that.mexp;
            pos = that.pos;
            sh1 = that.sh1;
            sh2 = that.sh2;
            id = that.id;
            mask = that.mask;
            for (int i = 0; i < 4; i++) {
                tbl[i] = that.tbl[i];
            }
            for (int i = 0; i < 16; i++) {
                p[i] = that.p[i];
            }
            for (int i = 0; i < 4; i++) {
                tmp_tbl[i] = that.tmp_tbl[i];
            }
            for (int i = 0; i < 16; i++) {
                tp[i] = that.tp[i];
            }
#if defined(USE_SHA)
            for (int i = 0; i < SHA_DIGEST_LENGTH * 2; i++) {
                sha1[i] = that.sha1[i];
            }
            sha1[0] = that.sha1[0];
#endif
        }

        void set_sha1(std::string& src) {
#if defined(USE_SHA)
            unsigned int sha1_length = SHA_DIGEST_LENGTH * 2 + 1;
            for (unsigned int i = 0; i < sha1_length && i < src.size(); i++) {
                sha1[i] = src[i];
            }
            sha1[sha1_length - 1] = 0;
#else
            MTToolBox::UNUSED_VARIABLE(&src);
            sha1[0] = 0;
#endif
        }

        const std::string getHeaderString() {
            return "sha1, mexp, type, id, pos, sh1, sh2, tbl_0,"
                " tbl_1, tbl_2, tbl_3,tmp_0, tmp_1, tmp_2, tmp_3,"
                " mask";
        }

        const std::string getParamString() {
            using namespace std;
            using namespace MTToolBox;
            stringstream ss;
            ss << '"' << sha1 << '"' << ",";
            ss << dec << mexp << ",";
            int length = 8;
            if (bit_size<T>() == 32) {
                ss << "uint32_t,";
                length = 8;
            } else {
                ss << "uint64_t,";
                length = 16;
            }
            ss << dec << id << ",";
            ss << dec << pos << ",";
            ss << dec << sh1 << ",";
            ss << dec << sh2 << ",";
            for (int i = 0; i < 4; i++) {
                ss << hex << setfill('0') << setw(8) << tbl[i] << ",";
            }
            for (int i = 0; i < 4; i++) {
                ss << hex << setfill('0') << setw(length) << tmp_tbl[i] << ",";
            }
            ss << hex << setfill('0') << setw(length) << mask << ",";
#if defined(DEBUG)
            for (int i = 0; i < 16; i++) {
                ss << hex << p[i] << ",";
            }
            for (int i = 0; i < 16; i++) {
                ss << hex << tp[i] << ",";
            }
#endif
            return ss.str();
        }

        bool equals(mtgp_param<T>& other) {
            if (mexp != other.mexp ||
                pos != other.pos ||
                sh1 != other.sh1 ||
                sh2 != other.sh2 ||
                id != other.id ||
                mask != other.mask) {
                return false;
            }
            for (int i = 0; i < 4; i++) {
                if (tbl[i] != other.tbl[i] ||
                    tmp_tbl[i] != other.tmp_tbl[i]) {
                    return false;
                }
            }
            return true;
        }

        int mexp;
        int pos;
        int sh1;
        int sh2;
        uint32_t id;
        T mask;
        T tbl[4];
        T p[16];
        T tmp_tbl[4];
        T tp[16];
#if defined(USE_SHA)
        char sha1[SHA_DIGEST_LENGTH * 2 + 1];
#else
        char sha1[1];
#endif
    };
}
#endif
