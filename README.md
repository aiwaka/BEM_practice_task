# BEM
境界要素法の練習用コード.
半径$1$の円周領域の境界上で与えられた関数$u(x)=x^3-3xy^2$に対して, 境界要素法を用いて共役な法線微分を計算し, また任意の内点での値を計算する.

# Compile command
LAPACKの導入については別途参照のこと. ```/usr/local/lib```にLAPACKがインストールされているなら
```
gfortran -g -O0 -o main main.f90 -I/usr/local/lib -llapack -lblas
```
と入力してコンパイルする.

# Execution
```
./main
```

# Input
1. 円周の分割数を入力.
2. 値を求める点を入力. 半角スペース区切りで2つ入力する. 例えば, ``` 0.25d0 0.5d0 ```など.