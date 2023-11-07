# 境界要素法でLaplace方程式を解くコード

## 概要

境界要素法の実装の練習用コード。

原点中心、半径1の円周領域の境界 $\Gamma$ を考え、Laplace方程式の境界値問題

$$
\Delta u = 0
$$

$$
u(x, y) = x^3-3xy^2 \quad \text{on} \quad \Gamma
$$

に対して境界要素法を用いて解を求める。
$u(x,y)$ に対して共役な法線微分を計算し、境界内部の領域の点での値を計算してその結果を出力する。

## コンパイル・実行

LAPACKの導入については別途参照のこと。
`set_env.sh.sample`を`set_env.sh`に改名して`source set_env.sh`を実行し、`make`を実行するとコンパイルされる。

```shell
./main.out
```

等で実行し、分割数を入力すると各点の問題の解の値が計算され、`plot.dat`に出力される。

## プロット

gnuplotを起動し`sp 'plot.dat' w p`などとすればプロットされる。
また、スクリプトを用いて`./plot.sh plot.gp`を実行してもよい。