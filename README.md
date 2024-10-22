#文件夹说明
+ dataSet
  + density:以密度筛选出的数据集
  + pc:通过堆积系数筛选出的数据集
  + D_chon.csv:chon类分子数据集
  + D_no2.csv:硝基类分子集合
  + data_run.py:数据筛选代码，包括测试集构建、分子指纹的最远点采样、数据整理等
  + testSet.csv:独立测试集
+ model 
  + generalizability_test:测试集与训练集差异性对模型泛化能力影响
  + model_density
  + model_D_chon
  + model_D_no2
  + model_D_pc  
  
  
文章模型使用 [chemProp](https://github.com/chemprop/chemprop) 得到，详细训练参数见model内的args.json文件

