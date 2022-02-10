'''
Decentralized scheme based on paper by Chotard, Dufour Sans, Gay, Phan and Pointcheval. 
https://eprint.iacr.org/2017/989.pdf
This scheme does not require a trusted party to generate keys. It is built on pairings and based on SXDH assumption.

Author: Liyao Xiong
Date: 2021.07.17 
'''

from contextlib import AsyncExitStack
import json
from test_mife_dynamic import parties_input_size
import numpy as np
import logging
import numpy as np
import hashlib
import math

import gmpy2 as gp
from crypto.dmcfe_utils import _random_with_seed

logger = logging.getLogger(__name__)


class DamgardDMCFEClient():
    def __init__(self, index=None):
        self.index = index  # 用户标签
        # if sec_param_config is not None and dlog is not None:
        #     self.p, self.q, self.r, self.g, self.sec_param = self._load_sec_param_config(
        #         sec_param_config)
        #     self.dlog_table = dlog['dlog_table']
        #     self.func_bound = dlog['func_bound']
        #     assert dlog['g'] == self.g, 'g in dlog table does not match g in sec param'
        # else:
        #     self.p, self.q, self.r = _param_generator(sec_param)
        #     self.g = _random_generator(sec_param, self.p, self.r)
        #     self.sec_param = sec_param
        #     self.dlog_table = None
        #     self.func_bound = None

    def load_sec_param_config(self, sec_param_config_file, parties, input_size):
        with open(sec_param_config_file, 'r') as infile:
            sec_param_dict = json.load(infile)

            p = gp.mpz(sec_param_dict['group']['p'])
            q = gp.mpz(sec_param_dict['group']['q'])
            h = gp.mpz(sec_param_dict['group']['h'])
            g = gp.mpz(sec_param_dict['g'])
            sec_param = sec_param_dict['sec_param']

        # 返回公共参数mpk
        return {'p': p,
                'q': q,
                'h': h,
                'g': g,
                'sec_param': sec_param,
                'parties': parties,
                'input_size': input_size}

    def setup(self, mpk):
        '''返回(ek, client_sec_key, client_pub_key)
        ek_i = s_i in Z_p^2：用户私钥的一部分
        client_sec_key in Z_p：用于后续生成用户私钥的另一部分T_i
        client_pub_key = g1 * client_sec_key in Z_p：用于后续生成用户私钥的另一部分T_i
        '''
        self.client_sec_key = _random_with_seed(gp.mpz(0), mpk['q'], 1)
        client_pub_key = gp.powmod(mpk['g'], self.client_sec_key, mpk['p'])

        return client_pub_key

    def set_share(self, client_pub_keys, mpk):
        ''' 生成T_i in Z_p^{2x2}，用户私钥的另一部分
            所有用户的T_i之和SUM(T_i) = 0
        '''
        input_size = mpk['input_size']
        parties = mpk['parties']
        total_num = parties * input_size
        # hash256 = hashlib.sha256()
        share = [[gp.mpz(0) for _ in range(input_size)]
                 for _ in range(parties)]

        for k in range(parties):
            if k == self.index:
                continue
            #clients[k]['client_pub_key'] = sec[k]*g1, clients[k]['client_sec_key']
                #sharedG1 = clients[k]['client_pub_key'] * clients[idx]['client_sec_key']
            #         = clients[idx]['client_sec_key'] * clients[k]['client_pub_key']
            sharedG1 = gp.powmod(
                client_pub_keys[k], self.client_sec_key, mpk['p'])
            # 生成sha256校验和字符串(len=64)，和sharedG1一一对应
            # 修改：直接用shareG1作为seed
            # hash256.update(str(sharedG1).encode('utf-8'))
            # checksum = int.from_bytes(hash256.digest(), 'big')
            
            add = _random_with_seed(
                gp.mpz(0), mpk['q'], total_num, int(sharedG1))

            if k < self.index:
                for i in range(total_num):
                    share[i//input_size][i % input_size] = share[i//input_size][i % input_size] + add[i]
            else:
                for i in range(total_num):
                    share[i//input_size][i % input_size] = share[i//input_size][i % input_size] - add[i]
            
        # debug
        # for i in range(total_num):
        #     logger.debug(share[i//input_size][i % input_size])

        self.share = share

    def gen_dam_sec_key(self, mpk):
        '''
        生成sec_key = (dam_sec_key, dam_pub_key, otp_key)
        i in range(input_len)
        dam_sec_key = (s, t), t[i]、s[i]为[2, p - 1)的随机数
        dam_pub_key[i] = ((g^s[i] mod p) * (h^t[i] mod p)) mod p
        otp_key[i] = [0, q)的随机数
        '''
        input_size = mpk['input_size']
        # s = _random_with_seed(gp.mpz(2), mpk['p'] - gp.mpz(1), input_size)
        # t = _random_with_seed(gp.mpz(2), mpk['p'] - gp.mpz(1), input_size)
        s = _random_with_seed(gp.mpz(2), mpk['q'], input_size)
        t = _random_with_seed(gp.mpz(2), mpk['q'], input_size)
        dam_pub_key = []
        otp_key = _random_with_seed(
            gp.mpz(0), mpk['q'], input_size)

        for i in range(input_size):
            y1 = gp.powmod(mpk['g'], s[i], mpk['p'])
            y2 = gp.powmod(mpk['h'], t[i], mpk['p'])
            dam_pub_key.append(gp.t_mod(gp.mul(y1, y2), mpk['p']))

        dam_sec_key = {'s': s, 't': t}
        self.sec_key = {'dam_sec_key': dam_sec_key,
                        'dam_pub_key': dam_pub_key, 'otp_key': otp_key}

    def encrypt(self, x, mpk):
        ''' 
        vector x_add_otp = (otp_key + x) mod q
        r = [1, p)随机数
        ct[0] = g^r mod p
        ct[1] = h^r mod p
        for i in range(input_size):
            t1 = dam_pub_key[i]^r mod p
            t2 = g^x_add_otp[i] mod p
            ct[i + 2] = (t1 * t2) mod p
        *即ct的元素个数为len+2
        '''
        input_size = mpk['input_size']
        #r = _random_with_seed(gp.mpz(1), mpk['p'], 1)
        r = _random_with_seed(gp.mpz(2), mpk['q'], 1)
        ct = []
        ct.append(gp.powmod(mpk['g'], r, mpk['p']))
        ct.append(gp.powmod(mpk['h'], r, mpk['p']))

        for i in range(input_size):
            x_add_otp = gp.t_mod(
                self.sec_key['otp_key'][i] + gp.mpz(x[i]), mpk['q'])
            t1 = gp.powmod(self.sec_key['dam_pub_key'][i], r, mpk['p'])
            t2 = gp.powmod(mpk['g'], x_add_otp, mpk['p'])
            ct.append(gp.t_mod(gp.mul(t1, t2), mpk['p']))

        return ct

    def derive_key_share(self, y, mpk):
        ''' DkeyGenShare(sk[i], vector y[idx], idx)
        fe_key_part = (otp_key_part, key_part)
        z1 = <otp_key, y[idx]>
        z2 = dot(share, y)
        otp_key_part = (z1 + z2) mod q

        key_part = (key1, key2)
        key1 = <dam_sec_key->s, y[idx]> mod (p - 1)
        key2 = <dam_sec_key->t, y[idx]> mod (p - 1)
        '''
        input_size = mpk['input_size']
        parties = mpk['parties']
        y_part = y[self.index]
        z1 = gp.mpz(0)
        key1 = gp.mpz(0)
        key2 = gp.mpz(0)
        for i in range(input_size):
            z1 = z1 + gp.mul(self.sec_key['otp_key'][i], gp.mpz(y_part[i]))
            key1 = key1 + \
                gp.mul(self.sec_key['dam_sec_key']['s'][i], gp.mpz(y_part[i]))
            key2 = key2 + \
                gp.mul(self.sec_key['dam_sec_key']['t'][i], gp.mpz(y_part[i]))
        # key_part = {'key1': gp.t_mod(key1, mpk['p'] - gp.mpz(1)),
        #             'key2': gp.t_mod(key2, mpk['p'] - gp.mpz(1))}
        key_part = {'key1': gp.t_mod(key1, mpk['q']),
                    'key2': gp.t_mod(key2, mpk['q'])}

        z2 = gp.mpz(0)
        for i in range(parties):
            for j in range(input_size):
                z2 = z2 + gp.mul(self.share[i][j], gp.mpz(y[i][j]))

        otp_key_part = gp.t_mod(z1 + z2, mpk['q'])
        return {'otp_key_part': otp_key_part, 'key_part': key_part}


class DamgardDMCFEServer():
    def __init__(self, dlog_table=None) -> None:
        self.dlog_table = dlog_table

    def decrypt(self, fe_key_parts, ciphers, mpk, y):
        ''' Dec(dk_y, l, c)
            z = SUM(fe_key_parts['otp_key_part']) mod q
            for i in range(parties):
                keys[i] = fe_key_parts[i]['key_part']

            alpha = 1
            for i in range(parties):
                num = 1
                for j in range(len):
                    t1 = ct[i][j+2]^key_part['key1'] mod p
                    num = (num * t1) mod p
                t1 = ct[i][0]^key_part['key1'] mod p
                t2 = ct[i][1]^key_part['key2'] mod p
                denom = t1 * t2
                denom_inv = invert(denom, p)
                alpha = (alpha * num * denom_inv) mod p

            z_exp = g^z mod p
            z_exp_inv = invert(z_exp, p)
            alpha = (alpha * z_exp_inv) mod p

            已知g^x mod p = alpha，解得alpha
        '''
        input_size = mpk['input_size']
        parties = mpk['parties']
        z = gp.mpz(0)
        key_part = []
        for i in range(parties):
            z = z + fe_key_parts[i]['otp_key_part']
            key_part.append(fe_key_parts[i]['key_part'])
        z = gp.t_mod(z, mpk['q'])

        alpha = gp.mpz(1)
        for i in range(parties):
            num = gp.mpz(1)
            for j in range(input_size):
                cy = gp.powmod(ciphers[i][j + 2],
                               gp.mpz(y[i][j]), mpk['p'])
                num = gp.t_mod(gp.mul(num, cy), mpk['p'])
            t1 = gp.powmod(ciphers[i][0], key_part[i]['key1'], mpk['p'])
            t2 = gp.powmod(ciphers[i][1], key_part[i]['key2'], mpk['p'])
            denom = gp.t_mod(t1 * t2, mpk['p'])
            denom_inv = gp.invert(denom, mpk['p'])
            alpha = gp.mul(alpha, denom_inv)
            alpha = gp.t_mod(gp.mul(alpha, num), mpk['p'])

        z_exp = gp.powmod(mpk['g'], z, mpk['p'])
        z_exp_inv = gp.invert(z_exp, mpk['p'])
        alpha = gp.t_mod(gp.mul(alpha, z_exp_inv), mpk['p'])

        return self._solve_dlog(alpha, mpk['g'], mpk['p'])

    def _solve_dlog(self, alpha, g, p):
        """
        Attempts to solve for the discrete log x, where g^x = h mod p via
        hash table.
        """
        if self.dlog_table is not None and alpha.digits() in self.dlog_table:
            return self.dlog_table[alpha.digits()]
        else:
            logger.debug(
                "did not find f in dlog table, may cost more time to compute")
            return self._solve_dlog_bsgs(g, alpha, p)

    def _solve_dlog_bsgs(self, g, h, p):
        """
        Attempts to solve for the discrete log x, where g^x = h mod p,
        via the Baby-Step Giant-Step algorithm.
        """
        m = math.ceil(math.sqrt(p-1))  # phi(p) is p-1, if p is prime
        # store hashmap of g^{1,2,...,m}(mod p)
        hash_table = {pow(g, i, p): i for i in range(m)}
        # precompute via Fermat's Little Theorem
        c = pow(g, m * (p-2), p)
        # search for an equivalence in the table. Giant Step.
        for j in range(m):
            y = (h * pow(c, j, p)) % p
            if y in hash_table:
                return j * m + hash_table[y]

        return None
