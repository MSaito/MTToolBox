double precision SIMD oriented Fast Mersenne Twister sample
===========================================================

dSFMTdc dSFMTのパラメータを生成する
usage:
./dSFMTdc [-s seed] [-v] [-c count] [-f outputfile] mexp

--verbose, -v        計算時間とか表示する。
--file, -f filename  探索されたパラメータはこのファイルに出力される。指定されなければ
	   	     標準出力に出る。
--count, -c count    出力件数
--seed, -s seed      乱数の種
--fixed, -x fixedSL1 sl1パラメータを指定された値に固定する。
mexp                 メルセンヌ指数

calc_equidist dSFMTのパラメータから均等分布次元を計算する
usage:
./calc_equidist [-v] mexp,pos1,sl1,msk1,msk2,fix1,fix2,parity1,parity2

--verbose, -v        各v毎のvビット精度均等分布次元を出力する。
外の引数はdSFMTdc の出力を使う。

test_linearity 線形性のテスト（デバッグ用）
./test_linearity
mexp,pos1,sl1,msk1,msk2,fix1,fix2,parity1,parity2
引数はdSFMTdc の出力を使う。

ちょっといい話かも知れない話
dSFMT のパラメータのシフト量は、instrincs関数で書くと、リテラルにしなければならない訳だが、
-x 引数によってシフト量を固定値にしてしまえば、シフト量が同じdSFMTのパラーメータを生成できる。
そうなればdSFMTの疑似乱数生成コードをコンパイルしなおさなくても、（シフト量以外の）パラメータ
の異なるdSFMTに対応できる。
ただし、固定したシフト量の値によっては、パラメータがまったく検索されないという可能性もある。
