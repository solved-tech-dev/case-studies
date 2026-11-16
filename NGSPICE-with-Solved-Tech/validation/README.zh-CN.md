# 原生方程采集验证

发布包构建的 ngspice 已完整运行五例。逐字节核对确认，原 CSR、每次 RHS 和 KLU 解向量，以及完整物理场数据均与原测试一致。物理场比较从工作点结果标题开始，独立于 raw 文件的日期与构建信息。详见 [capture-validation.json](capture-validation.json)。
