import logging
from crypto.damgard_dmcfe import DamgardDMCFEServer, DamgardDMCFEClient
from crypto.dmcfe_utils import generate_config_files, load_dlog_table_config

import random
import os
import numpy as np
from nn.utils import timer
import gmpy2 as gp

logger = logging.getLogger(__name__)

sec_param_config = 'config/sec_param.json'
dlog_table_config = 'config/dlog.json'


def test_crypto_config():
    func_value_bound = 100000000
    #func_value_bound = 10000000
    sec_param = 256
    max_len = 784

    if not (os.path.exists(sec_param_config) and os.path.exists(dlog_table_config)):
        logger.debug(
            "could not find the crypto config file, generate a new one")
        generate_config_files(sec_param, sec_param_config,
                              dlog_table_config, func_value_bound, max_len)
    else:
        logger.info("the crypto config file already existed")
    return sec_param_config, dlog_table_config


def test_damgard_dmcfe_with_config():
    logger.info("testing the correctness of damgard dmcfe using config file.")

    max_test_value = 256 # 最大值
    input_len = 200 #每个用户向量长度
    parties = 5 # 用户数
    x_mat = []
    x_vec = []
    y_mat = []
    y_vec = []
    for idx in range(parties):
        x_mat.append([random.randint(0, max_test_value)
                     for _ in range(input_len)])
        y_mat.append([random.randint(0, max_test_value)
                     for _ in range(input_len)])
        x_vec = x_vec + x_mat[idx]
        y_vec = y_vec + y_mat[idx]

    logger.info('loading dlog configuration ...')
    with timer('load dlog config, cost time:', logger) as t:
        dlog = load_dlog_table_config(dlog_table_config)
    logger.info('load dlog configuration DONE')
    logger.debug("func_bound = %s" % dlog['func_bound'])

    with timer('Setup phase, cost time:', logger) as t:
        mpk = DamgardDMCFEClient().load_sec_param_config(sec_param_config, parties, input_len)

        clients = []
        client_pub_keys = []
        for i in range(parties):
            client = DamgardDMCFEClient(i)
            clients.append(client)
            client_pub_keys.append(client.setup(mpk))
            client.gen_dam_sec_key(mpk)

        fe_key_part = []
        ct = []
        for i in range(parties):
            clients[i].set_share(client_pub_keys, mpk)

    with timer('DKeyShare phase, cost time:', logger) as t:
        for i in range(parties):
            fe_key_part.append(clients[i].derive_key_share(y_mat, mpk))

    with timer('Encrypt phase, cost time:', logger) as t:
        for i in range(parties):
            ct.append(clients[i].encrypt(x_mat[i], mpk))

    # 测试share是否相加为0
    # sum_mat = [[gp.mpz(0) for _ in range(input_len)] for _ in range(parties)]
    # for i in range(parties):
    #     for j in range(parties):
    #         for k in range(input_len):
    #             sum_mat[j][k] = sum_mat[j][k] + clients[i].share[j][k]
    # for j in range(parties):
    #     for k in range(input_len):
    #         logger.debug("sum_mat = %s\t", sum_mat[j][k])
    
    logger.debug("x: {}".format(x_vec))
    logger.debug("y: {}".format(y_vec))
    logger.debug('original dot product <x,y>: {}'.format(
        sum(np.array(x_vec) * np.array(y_vec))))

    with timer('DKeyComb and Decrypt phase, cost time:', logger) as t:
        server = DamgardDMCFEServer(dlog['dlog_table'])
        xy_cal = server.decrypt(fe_key_part, ct, mpk, y_mat)

    logger.debug('decrypted dot product <x,y>: {}'.format(xy_cal))
    assert xy_cal == sum(np.array(x_vec) * np.array(y_vec))
