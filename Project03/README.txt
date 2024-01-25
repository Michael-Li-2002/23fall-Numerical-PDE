Project03 李佳 2100010793

主程序: project03_main.m:            一般情形, 时间步长固定
            project03_specialcase.m:  特殊情形(\theta=1/2,1/2-1/(12\mu), 使用data_2,data_3(非光滑初值)), 
                                                      先用小时间步长, 再用大步长.
函数   : theta_HeatEq:                   求解热方程模型问题
            tridiagsolver:                     求解三对角矩阵线性方程, 并给出运算量
            Linferr:                               计算某时间层的L^\infty误差
            L2err:                                 计算某时间层的L^2误差
            data_1:                              光滑初值及真解
            data_2:                              连续但导数间断初值及真解(有截断)
            data_3:                              分片连续初值及真解(有截断)

每次运行对所有N=8,16,32,64,128计算, 
需要在主程序 Parameter 中给定 :
一个终止时间t, 每个N对应的\tau, 一个\theta, （若是project03_specialcase.m, 则需给出分界时间t_0前后各自的\tau,\theta）
运行后输出N, 运算量, L^\infty, L^2 误差的表格, 
并作出 L^2 误差分别与 N 、运算量的双对数图