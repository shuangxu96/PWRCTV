import os
import numpy as np
from scipy.io import loadmat
from PIL import Image


def make_noisy(dataset, casenum):
    """
    生成带噪 HSI 数据，对应 MATLAB 函数 make_noisy。

    参数:
        dataset : str  - 数据集名称 ('Florence' 或 'Milan')
        casenum : int  - 1~5，噪声类型

    返回:
        Nhsi : ndarray (256, 256, B) 带噪 HSI
        Ohsi : ndarray (256, 256, B) 干净 HSI (归一化)
        Pan  : ndarray (256, 256)    全色图像 (归一化)
    """
    root = '.'

    # (1) 读取干净的 HSI
    mat_path = os.path.join(root, 'data', dataset, f'{dataset}.mat')
    S = loadmat(mat_path)

    hsi = S['hsi']
    Ohsi = hsi.astype(np.float64) / 65535.0

    # (2) 读取 PAN 图像
    if 'pan' in S and S['pan'] is not None and S['pan'].size > 0:
        Pan = S['pan'].astype(np.float64) / 65535.0
        if Pan.ndim == 3 and Pan.shape[2] == 1:
            Pan = Pan[:, :, 0]
    else:
        # 如果没有 pan 图像，用波段平均模拟
        Pan = np.mean(Ohsi, axis=2)

    M, N, B = Ohsi.shape

    # 重采样到 256x256
    if M != 256 or N != 256:
        Ohsi = _imresize(Ohsi, (256, 256))
        Pan = _imresize(Pan, (256, 256))

    M, N, B = Ohsi.shape

    # (3) 设置随机种子
    ds_names = ['florence', 'milan']
    ds_idx = -1
    for k, name in enumerate(ds_names):
        if dataset.lower() == name:
            ds_idx = k + 1  # MATLAB find 是 1-based
            break
    if ds_idx == -1:
        ds_idx = 0
    np.random.seed(1000 * ds_idx + casenum)

    # (4) 噪声
    Nhsi = Ohsi.copy()

    if casenum == 1:
        sig = 10.0 / 255.0                                  # iid
    elif casenum in (2, 3, 4, 5):
        sig = (np.random.rand() * 25 + 5) / 255.0           # non-iid
    else:
        sig = 0.0

    # (4a) 高斯噪声
    for i in range(B):
        Nhsi[:, :, i] = Nhsi[:, :, i] + np.random.randn(M, N) * sig

    # (4b) case 3&5：脉冲噪声
    if casenum in (3, 5):
        ipb = np.random.permutation(B)[:int(np.ceil(0.333 * B))]
        for i in ipb:
            p = np.random.rand() * 0.25 + 0.05              # 5%~30%
            Nhsi[:, :, i] = _salt_pepper(Nhsi[:, :, i], p)

    # (4c) case 4&5：条纹噪声
    if casenum in (4, 5):
        stb = np.random.permutation(B)[:int(np.ceil(0.333 * B))]
        for i in stb:
            s = np.random.rand() * 0.25 + 0.05
            linenum = int(np.ceil(N * s))
            lineloc = np.random.permutation(N)[:linenum]    # 不重复的列位置
            # MATLAB: ceil(N*rand(1,linenum)) 可能重复；这里保持简单随机（1-based -> 0-based）
            # 若需严格一致，可改为 np.random.randint(0, N, linenum)
            t = np.random.rand(len(lineloc)) * 0.5 - 0.25
            Nhsi[:, lineloc, i] = Nhsi[:, lineloc, i] - t

    # (5) 裁剪到 [0,1]
    np.clip(Nhsi, 0.0, 1.0, out=Nhsi)

    return Nhsi, Ohsi, Pan


def _imresize(img, size):
    """使用 PIL 对 float 图像做双线性插值缩放（保持 [0,1] 范围）。"""
    h, w = size
    if img.ndim == 2:
        img_min = img.min()
        img_max = img.max()
        rng = img_max - img_min
        if rng == 0:
            rng = 1.0
        arr = (img - img_min) / rng
        arr = np.clip(arr, 0, 1)
        pil = Image.fromarray((arr * 255).astype(np.uint8), mode='L')
        pil = pil.resize((w, h), Image.BILINEAR)
        out = np.asarray(pil, dtype=np.float64) / 255.0
        return out * rng + img_min
    else:
        B = img.shape[2]
        out = np.zeros((h, w, B), dtype=np.float64)
        for i in range(B):
            out[:, :, i] = _imresize(img[:, :, i], size)
        return out


def _salt_pepper(img, p):
    """与 MATLAB imnoise(..., 'salt & pepper', p) 一致。"""
    out = img.copy()
    r = np.random.rand(*img.shape)
    salt = r < (p / 2.0)
    pepper = (r >= (p / 2.0)) & (r < p)   # 用 & 做逐元素布尔运算
    out[salt] = 1.0
    out[pepper] = 0.0
    return out