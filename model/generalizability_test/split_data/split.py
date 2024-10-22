
from rdkit.Chem.PandasTools import pd
from sklearn.metrics import r2_score, mean_absolute_error

from chemprop.train import make_predictions

from chemprop.args import PredictArgs


def splitData():
    reader = open(r"./35250.csv", "r")
    base = reader.readlines()
    title = base[0]
    for i in range(int(len(base) / 1000)):
        writer = open("./%s.csv" % str(i + 1), "x")
        writer.write(title)
        writer.writelines(base[i * 1000 + 1:(i + 1) * 1000 + 1])
        writer.flush()
        writer.close()


# def testDataValue():
#     for i in range(1, 36):
#         for j in range(1, 36):
#             predict_args = ['--checkpoint_dir', "D:\pycharmProject\chemprop-master\model\split/%s" % str(i)
#                 , '--test_path', "D:\pycharmProject\chemprop-master/model/split/%s" % str(j) + "/fold_0/test_full.csv"
#                 , '--preds_path', "D:\pycharmProject\chemprop-master/predict/split/%s_%s.csv"%(str(i),str(j))
#                 , "--smiles_columns", "smiles"
#                             ]
#             prediction_args = PredictArgs().parse_args(predict_args)
#             make_predictions(args=prediction_args)


def scores():
    title = [str(i) for i in range(0,37)]
    writer1 = open("./score_r2_2.csv","x")
    writer1.write(",".join(title)+"\n")
    writer2 = open("./score_mae_2.csv", "x")
    writer2.write(",".join(title)+"\n")
    for i in range(1, 36):
        density_path = "D:\pycharmProject\chemprop-master/model/split/%s" % str(i) + "/fold_0/test_full.csv"
        density = pd.read_csv(density_path, delim_whitespace=False, header=0)["Density"]
        score_r2 = [str(i)]
        score_mae = [str(i)]
        for j in range(1, 37):
            pred_path = "D:\pycharmProject\chemprop-master/predict/split/%s_%s.csv"%(str(j),str(i))
            pred = pd.read_csv(pred_path, delim_whitespace=False, header=0)["Density"]
            r2 = r2_score(density, pred)
            mae = mean_absolute_error(density, pred)
            score_r2.append(str(r2))
            score_mae.append(str(mae))
        print(score_r2)
        print(score_mae)
        writer1.write(",".join(score_r2)+"\n")
        writer2.write(",".join(score_mae)+"\n")


def scoreOverAll():
    for i in range(36, 37):
        for j in range(1, 36):
            predict_args = ['--checkpoint_dir', r"D:\pycharmProject\chemprop-master\model\train_35250"
                , '--test_path', "D:\pycharmProject\chemprop-master/model/split/%s" % str(j) + "/fold_0/test_full.csv"
                , '--preds_path', "D:\pycharmProject\chemprop-master/predict/split/36_%s.csv" % (str(j))
                , "--smiles_columns", "smiles"
                            ]
            prediction_args = PredictArgs().parse_args(predict_args)
            make_predictions(args=prediction_args)

if __name__ == '__main__':
    # splitData()
    # testDataValue()
    scores()
    # scoreOverAll()
    # print("s")