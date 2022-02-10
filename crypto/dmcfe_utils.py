import random
import zlib
import base64
import json
import logging
from datetime import datetime

import gmpy2 as gp

logger = logging.getLogger(__name__)

def _random_with_seed(minimum, maximum, num, seed=None):
    '''
    相同种子生成的随机序列相同
    生成[min, max)随机数
    '''
    if seed == None:
        state = gp.random_state(int('{0:%Y%m%d%H%M%S%f}'.format(datetime.now())))
    else:
        # 不指定seed时，以当前时间为seed，保证每次随机的数字不一样
        state = gp.random_state(seed)
    ret = []

    for _ in range(num):
        # mpz_random：根据state生成[0, maximum)的随机数
        r = gp.mpz_random(state, maximum)
        while r < minimum:
            r = gp.mpz_random(state, maximum)
        ret.append(r)

    if num == 1:
        return ret[0]
    else:
        return ret

def _param_generator(bits):
    '''
    生成参数
    '''
    rand_function = random.SystemRandom()
    while True:
        q = gp.mpz(rand_function.getrandbits(bits))
        q = gp.next_prime(q)
        if gp.is_prime(q) == False:
            continue
        p = gp.mul(q, gp.mpz(2)) + gp.mpz(1)
        if gp.is_prime(p) == True and p.bit_length() == bits:
            break
    
    logger.info("Generate prime q and p success")

    while True:
        r = _random_with_seed(gp.mpz(3), p, 1)
        g = gp.powmod(r, gp.mpz(2), p)
        g_inv = gp.invert(g, p)
        # additional checks to avoid some known attacks
        if gp.t_mod(p - 1, g) != 0 and gp.t_mod(p-1, g_inv) != 0:
            logger.info("Generate g success")
            break

    return g, p, q


def generate_config_files(sec_param, sec_param_config, dlog_table_config, func_bound, max_len):
    '''
    g: 循环群Z_P的生成元：g^q = 1 (mod p)
    h: 循环群Z_P的生成元：h^q = 1 (mod p)
    p: 模数, p = 2 * q + 1
    q: g和h的乘法阶
    '''
    g, p, q = _param_generator(sec_param)
    if gp.mpz(func_bound ** 2 * max_len) > q:
        logger.debug("bound error")
        return

    
    while True:
        # h generated in the following way is always a generator with order q
        # r = _random_with_seed(gp.mpz(1), p, sec_param, 1)
        r = _random_with_seed(gp.mpz(2), q, 1)
        h = gp.powmod(g, r, p)
        # additional checks to avoid some known attacks
        h_inv = gp.invert(g, p)
        if gp.t_mod(p - gp.mpz(1), h) != 0 and gp.t_mod(p - gp.mpz(1), h_inv) != 0:
            logger.info("Generate h success")
            break

    group_info = {
        'p': gp.digits(p), 
        'q': gp.digits(q), 
        'h': gp.digits(h)
    }
    sec_param_dict = {'g': gp.digits(
        g), 'sec_param': sec_param, 'group': group_info}

    with open(sec_param_config, 'w') as outfile:
        json.dump(sec_param_dict, outfile)
    logger.info(
        'Generate secure parameters config file successfully, see file %s' % sec_param_config)

    dlog_table = dict()
    bound = func_bound + 1
    for i in range(bound):
        dlog_table[gp.digits(gp.powmod(g, i, p))] = i
    for i in range(-1, -bound, -1):
        dlog_table[gp.digits(gp.powmod(g, i, p))] = i

    dlog_table_dict = {
        'g': gp.digits(g),
        'func_bound': func_bound,
        'dlog_table': dlog_table
    }

    with open(dlog_table_config, 'w') as outfile:
        # outfile.write(_json_zip(dlog_table_dict))
        json.dump(dlog_table_dict, outfile)
    logger.info(
        'Generate dlog table config file successfully, see file %s' % dlog_table_config)


def load_sec_param_config(sec_param_config_file):
    with open(sec_param_config_file, 'r') as infile:
        sec_param_dict = json.load(infile)

        p = gp.mpz(sec_param_dict['group']['p'])
        q = gp.mpz(sec_param_dict['group']['q'])
        h = gp.mpz(sec_param_dict['group']['h'])
        g = gp.mpz(sec_param_dict['g'])
        sec_param = sec_param_dict['sec_param']

    return p, q, h, g, sec_param


def load_dlog_table_config(dlog_table_config_file):
    with open(dlog_table_config_file, 'r') as infile:
        # config_content = infile.read()
        # store_dict = _json_unzip(config_content)
        store_dict = json.load(infile)

        dlog_table = store_dict['dlog_table']
        func_bound = store_dict['func_bound']
        g = gp.mpz(store_dict['g'])

    return {
        'dlog_table': dlog_table,
        'func_bound': func_bound,
        'g': g
    }
