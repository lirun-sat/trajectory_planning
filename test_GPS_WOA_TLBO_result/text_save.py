def text_save(filename, data_1, data_2, data_3, data_4):  # filename 为写入文件的路径(用字符串表示)，data 为要写入的数据列表.
    file = open(filename, 'a')  # a 表示在原来的基础上增加要写入的内容

    for i in range(len(data_1)):
        s_1 = str(data_1[i]).replace('[', '').replace(']', '')  # 去除[],这两行按数据不同，可以选择
        s_1 = s_1.replace("'", '').replace(',', '') + '\n'   # 去除单引号，逗号，每行末尾追加换行符
        file.write(s_1)
    file.write("\n")

    for i in range(len(data_2)):
        s_2 = str(data_2[i]).replace('[', '').replace(']', '')  # 去除[],这两行按数据不同，可以选择
        s_2 = s_2.replace("'", '').replace(',', '') + '\n'   # 去除单引号，逗号，每行末尾追加换行符
        file.write(s_2)
    file.write("\n")
    for i in range(len(data_3)):
        s_3 = str(data_3[i]).replace('[', '').replace(']', '')  # 去除[],这两行按数据不同，可以选择
        s_3 = s_3.replace("'", '').replace(',', '') + '\n'  # 去除单引号，逗号，每行末尾追加换行符
        file.write(s_3)
    file.write("\n")
    for i in range(len(data_4)):
        s_4 = str(data_4[i]).replace('[', '').replace(']', '')  # 去除[],这两行按数据不同，可以选择
        s_4 = s_4.replace("'", '').replace(',', '') + '\n'  # 去除单引号，逗号，每行末尾追加换行符
        file.write(s_4)
    file.write("\n")

    file.close()
    print("保存文件成功")
