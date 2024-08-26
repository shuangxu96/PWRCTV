# Hyperspectral Pan-denoising by PWRCTV

[Shuang Xu](https://shuangxu96.github.io/), [Qiao Ke](https://teacher.nwpu.edu.cn/qiaoke.html), [Jiangjun Peng](https://teacher.nwpu.edu.cn/pengjj), [Xiangyong Cao](https://gr.xjtu.edu.cn/en/web/caoxiangyong), [Zixiang Zhao](https://zhaozixiang1228.github.io/)

Northwestern Polytechnical University, and Xi'an Jiaotong University

-------------------------------------------
[![arXiv](https://img.shields.io/badge/arXiv-Paper-<COLOR>.svg)](https://arxiv.org/pdf/2407.06064)
[![Free](https://img.shields.io/badge/free_for_non_commercial_use-brightgreen)](#-license)

This repository contains the MATLAB code for the paper "Pan-denoising: Guided Hyperspectral Image Denoising via Weighted Represent Coefficient Total Variation". The code implements the proposed PWRCTV method for hyperspectral image denoising using panchromatic image guidance.

## Pan-denoising
With the recent launch of satellites equipped with both hyperspectral and PAN sensors, such as PRISMA (PRecursore IperSpettrale della Missione Applicativa) and XG3 (XIGUANG-003), a new opportunity has emerged. PAN images, due to their imaging mechanism, are less noisy than HSI but still exhibit similar textures. As depicted in following figure (b), this paper therefore aims to investigate PAN image-guided HSI denoising, which is referred to as _pan-denoising_. This problem arises from two primary aspects:
- Despite the significant advancements in hyperspectral imaging techniques, the HSIs captured by recent satellite sensors still suffer from noticeable noise. _Pan-denoising_ presents an important and novel approach to enhance HSI quality.
- Substantial research has been conducted on hyperspectral pan-sharpening, which assumes that HSIs are noise-free. However, this assumption does not hold in practice. Following pan-sharpening, a denoising step is still required. _Pan-denoising_ would lead to a more robust image preprocessing result.

Compared with the traditional HSI denoising paradigm, _pan-denoising_ incorporates an additional regularization term derived from external prior knowledge:

$$\min_{\mathcal{X}}\, \mathscr{L}\left( \mathcal{X},\mathcal{Y} \right)  + \lambda \mathscr{R}\left( \mathcal{X} \right)  +\tau \mathscr{E}\left( \mathcal{X},\mathbf{P} \right),$$

where $\mathscr{E}\left( \mathcal{X},\mathbf{P} \right)$ characterizes the external prior knowledge, $\mathbf{P}\in\mathbb{R}^{M\times N}$ is the PAN image, and $\tau$ controls the regularization strength. Nevertheless, designing an appropriate regularization term to effectively utilize the complementary information from PAN images remains a significant challenge.


<div align=center><img  src="pandenoising.svg"/></div>

## Contents
* `PWRCTV.m`: The implementation of the PWRCTV denoising algorithm.
* `Table2_3_Fig5_6_PRISMA.m`: Functions for evaluating the denoising performance on Florence and Milan datasets to reproduce Tables 2-3 and Figs. 5-6 of the manuscript.
* `Fig7_8_XG3.m`: Functions for evaluating the denoising performance on Beijing and Yulin datasets to reproduce Figs. 7-8 of the manuscript.
* `Table_4_step1.m`: Functions for denoising on the Urban dataset.
* `Table_4_step2.py`: Functions for classification on the denoised Urban dataset.
* `datasets/`: Dataset download links.
  
## Usage
1. **Data Download**: Download the datasets from the provided links and put them in the  `data/` folder.
2. **Experiments on simulated datasets**: Run `Table2_3_Fig5_6_PRISMA.m`.
3. **Experiments on real-world datasets**: Run `Fig7_8_XG3.m`.
4. **Visualize results**: The denoised images and corresponding metrics will be saved in the `results/` folder.

## Dependencies
* MATLAB

## Citation
If you use this code in your research, please cite the corresponding paper:
```
@article{PWRCTV,
  author       = {Shuang Xu and
                  Qiao Ke and
                  Jiangjun Peng and
                  Xiangyong Cao and
                  Zixiang Zhao},
  title        = {Pan-denoising: Guided Hyperspectral Image Denoising via Weighted Represent Coefficient Total Variation},
  journal      = {{IEEE} Trans. Geosci. Remote. Sens.},
  volume       = {62},
  pages        = {1--14},
  year         = {2024},
  doi          = {https://doi.org/10.1109/TGRS.2024.3450888},
}
```

## Contact
If you have any questions or need further assistance, please contact Shuang Xu at xs@nwpu.edu.cn
