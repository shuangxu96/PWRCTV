# Hyperspectral Pan-denoising by PWRCTV
This repository contains the MATLAB code for the paper "Pan-denoising: Guided Hyperspectral Image Denoising via Weighted Represent Coefficient Total Variation". The code implements the proposed PWRCTV method for hyperspectral image denoising using panchromatic image guidance.

## Pan-denoising
With the recent launch of satellites equipped with both hyperspectral and PAN sensors, such as PRISMA (PRecursore IperSpettrale della Missione Applicativa) and XG3 (XIGUANG-003), a new opportunity has emerged. PAN images, due to their imaging mechanism, are less noisy than HSI but still exhibit similar textures. As depicted in Fig. \ref{fig:toc}(b), this paper therefore aims to investigate PAN image-guided HSI denoising, which is referred to as \textit{pan-denoising}. This problem arises from two primary aspects:
	
		\item Despite the significant advancements in hyperspectral imaging techniques, the HSIs captured by recent satellite sensors still suffer from noticeable noise. \textit{Pan-denoising} presents an important and novel approach to enhance HSI quality.
		\item Substantial research has been conducted on hyperspectral pan-sharpening, which assumes that HSIs are noise-free. However, this assumption does not hold in practice. Following pan-sharpening, a denoising step is still required. \textit{Pan-denoising} would lead to a more robust image preprocessing result.
	\end{enumerate}
	
	Compared with the traditional HSI denoising paradigm as depicted in Eq. \eqref{eq:HSID}, pan-denoising incorporates an additional regularization term derived from external prior knowledge:
	$$
		\min_{\mathcal{X}}\, \mathscr{L}\left( \mathcal{X},\mathcal{Y} \right)  + \lambda \mathscr{R}\left( \mathcal{X} \right)  +\tau \mathscr{E}\left( \mathcal{X},\mathbf{P} \right),
	$$
	where $\mathscr{E}\left( \mathcal{X},\mathbf{P} \right)$ characterizes the external prior knowledge, $\mathbf{P}\in\mathbb{R}^{M\times N}$ is the PAN image, and $\tau$ controls the regularization strength. Nevertheless, designing an appropriate regularization term to effectively utilize the complementary information from PAN images remains a significant challenge.
![Pan-denoising.](pandenoising.svg)

## Contents
* `main.m`: The main function for running the denoising experiment.
* `PWRCTV.m`: The implementation of the PWRCTV denoising algorithm.
* `data_generator.m`: Functions for generating synthetic noisy hyperspectral images.
* `evaluate.m`: Functions for evaluating the denoising performance using various metrics.
* `datasets/`: Folder containing the synthetic and real-world datasets used in the experiments.
  
## Usage
1. **Setup**: Download the code and datasets from the repository.
2. **Run the main function**: Execute `main.m` to run the denoising experiment on the synthetic and real-world datasets. You can modify the parameters in the main function to suit your needs.
3. **Visualize results**: The denoised images and corresponding metrics will be saved in the `results/` folder.

## Dependencies
* MATLAB


