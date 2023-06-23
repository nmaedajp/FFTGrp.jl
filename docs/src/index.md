```@meta
CurrentModule = FFTGrp
```

# FFTGrp

## FFTを行う際の周辺関数

### ex_trend(y)
* データ y のリニアトレンドを取り除く

### dws(L; w_type="ct", base=1)
* L個のデータに対する，データウィンドウを返す．
* 初期設定は，cosine taper．
* データウィンドウを計算するとき，データの順序$j$を$-1 \leq t \leq 1$ の$t$に変換しているが，base = 1 は，配列の最初が１から始まるという意味で，$j=1$を$t=-1$に変換している．
* ウィンドウのタイプ
  * w_type = "p" || "parzen"    : parzen ウィンドウ
  * w_type = "t" || "triangle"  : 三角形ウィンドウ
  * w_type = "r" || "rectangle" : 長方形ウィンドウ
  * w_type = "b" || "boxcar"    : 長方形ウィンドウ
  * w_type = "w" || "welch"     : Welch ウィンドウ
  * w_type = "ct" || "cosine taper" : cosine taper
  * w_type = "han" || "hanning"  : ハニングウィンドウ

### swin_sq(Xf, band, N, dt)
* フーリエ変換に対するスペクトルウィンドウ
  * Xf：フーリエ変換
  * band: バンド幅
  * N：Xf のデータ数
  * dt：サンプリング周期
### swin(Pf, band, N, dt)
* パワースペクトルに対するスペクトルウィンドウ
  * Pf：パワースペクトル（two sided）
  * band: バンド幅
  * N：Xf のデータ数
  * dt：サンプリング周期

Documentation for [FFTGrp](https://github.com/nmaedajp/FFTGrp.jl).

```@index
```

```@autodocs
Modules = [FFTGrp]
```
