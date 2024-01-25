Project04 李佳 2100010793

主程序: project04.m       
函数   : Upwind:                             求解双曲方程的非守恒迎风格式
            Upwind_conservative:       求解双曲方程的守恒迎风格式
            LW_conservative:              求解双曲方程的 Lax-Wendroof 格式
            Linferr:                               计算某时间层的 L^\infty 误差
            L2err:                                 计算某时间层的 L^2 误差
            Burgers:                             Burgers方程通量函数、初值及真解(弱解)

在设置中更改: 
    option.h( x 方向步长 ), option.t( 终止时间 t ), option.tau( 时间步长)
    option.FDS( 有限差分格式: Upwind, Upwind_conservative, LW_conservative )
可输出对应时刻的数值解图像及其与弱解的 L^\infty, L^2 误差( 但由于解间断, L^\infty 误差的意义不大 )
注释部分用于批量生成上机报告中的数据