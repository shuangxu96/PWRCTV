# -*-coding=utf-8 -*-
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
import numpy as np
import scipy.io as sio
import random
import matplotlib.pyplot as plt
import os
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier

# label
label = sio.loadmat(r'data\Urban\end4_groundTruth.mat')
label = label['A']# pair with your dataset
label = np.reshape(label, [-1,307,307])
label = np.transpose(label, [2,1,0])
label = label.argmax(-1)+1
sample_num = 200

oa_ = []
aa_ = []
kappa_ = []
ac_ = []
# test methods 
method = ['Noisy', 'TCTV', 'LMHTV',  'LTHTV', 'LRTV', 'NGMeet', 'RCTV', 'WNLRATV', 
          'BALMF', 'CTV', 'RCILD', 'HLRTF', 'PWRCTV']
for index in range(len(method)):
    im = sio.loadmat(r'result\urban\%s.mat'%(method[index]))
    if method[index]=='Noisy':
        im = im['Nhsi']
    else:
        im = im['output']
        im = np.float64(im)/255


    # normalization
    im = (im - float(np.min(im)))
    im = im/np.max(im)

    # prepare dataset
    deepth = im.shape[2]
    classes = np.max(label)
    test_bg = False

    data_pos = {}
    train_pos = {}
    test_pos = {}
    
    for i in range(1,classes+1):
        data_pos[i]=[]
        train_pos[i]=[]
        test_pos[i] = []
    
    for i in range(label.shape[0]):
        for j in range(label.shape[1]):
            for k in range(1,classes+1):
                if label[i,j]==k:
                    data_pos[k].append([i,j])
                    continue
            
    random.seed(0)
    for i in range(1,classes+1): 
        indexies = random.sample(range(len(data_pos[i])),sample_num)
        for k in range(len(data_pos[i])):
            if k not in indexies:
                test_pos[i].append(data_pos[i][k])
            else:
                train_pos[i].append(data_pos[i][k])
    
    train = []
    train_label = []
    test = []
    test_label = []

    for i in range(1,len(train_pos)+1):
        for j in range(len(train_pos[i])):
            row,col = train_pos[i][j]
            train.append(im[row,col])
            train_label.append(i)
    
    for i in range(1,len(test_pos)+1):
        for j in range(len(test_pos[i])):
            row,col = test_pos[i][j]
            test.append(im[row,col])
            test_label.append(i)

    if not os.path.exists(r'result\urban\classification'):
        os.makedirs(r'result\urban\classification')

    # clf = SVC(C=100,kernel='rbf',gamma=1)
    clf = RandomForestClassifier(
        max_depth=5, n_estimators=10, max_features=1, random_state=42
    )
    train = np.asarray(train)
    train_label = np.asarray(train_label)
    clf.fit(train,train_label)
    C = np.max(label)
    
    # metrics
    matrix = np.zeros((C,C))
    for i in range(len(test)):
        r = clf.predict(test[i].reshape(-1,len(test[i])))
        matrix[r-1,test_label[i]-1] += 1
    
    ac_list = []
    for i in range(len(matrix)):
        ac = matrix[i, i] / sum(matrix[:, i])
        ac_list.append(ac)
        print(i+1,'class:','(', matrix[i, i], '/', sum(matrix[:, i]), ')', ac)
    print('confusion matrix:')
    print(np.int_(matrix))
    print('total right num:', np.sum(np.trace(matrix)))
    print('total test num:',np.sum(matrix))
    accuracy = np.sum(np.trace(matrix)) / np.sum(matrix)
    print('Overall accuracy:', accuracy)
    oa_.append(accuracy)
    # kappa
    kk = 0
    for i in range(matrix.shape[0]):
        kk += np.sum(matrix[i]) * np.sum(matrix[:, i])
    pe = kk / (np.sum(matrix) * np.sum(matrix))
    pa = np.trace(matrix) / np.sum(matrix)
    kappa = (pa - pe) / (1 - pe)
    ac_list = np.asarray(ac_list)
    aa = np.mean(ac_list)
    print('Average accuracy:',aa)
    print('Kappa:', kappa)
    aa_.append(aa)
    kappa_.append(kappa)
    ac_.append(ac_list)
    
    sio.savemat(os.path.join(r'result\urban\classification', '%s.mat'%(method[index])), 
                {'oa': accuracy,'aa':aa,'kappa':kappa,'ac_list':ac_list,'matrix':matrix})
    
    # classification map
    iG = np.zeros((label.shape[0],label.shape[1]))
    for i in range(label.shape[0]):
        for j in range(label.shape[1]):
            if label[i,j] == 0:
                if test_bg:
                    iG[i,j] = (clf.predict(im[i,j].reshape(-1,len(im[i,j]))))
                else:
                    iG[i,j]=0
            else:
                iG[i,j] = (clf.predict(im[i,j].reshape(-1,len(im[i,j]))))
    if test_bg:
        iG[0,0] = 0


    de_map = iG[::-1]
    fig, _ = plt.subplots()
    height, width = de_map.shape
    fig.set_size_inches(width/100.0, height/100.0)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.subplots_adjust(top=1,bottom=0,left=0,right=1,hspace=0,wspace=0)
    plt.axis('off')
    plt.axis('equal')
    plt.pcolor(de_map, cmap='jet')
    plt.savefig(os.path.join(r'result\urban\classification', '%s_map.png'%(method[index])),format='png',dpi=600)#bbox_inches='tight',pad_inches=0)
    plt.close()
    print('decode map get finished')

ac_ = np.array(ac_)

de_map = label[::-1]
fig, _ = plt.subplots()
height, width = de_map.shape
fig.set_size_inches(width/100.0, height/100.0)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.subplots_adjust(top=1,bottom=0,left=0,right=1,hspace=0,wspace=0)
plt.axis('off')
plt.axis('equal')
plt.pcolor(de_map, cmap='jet')
plt.savefig(os.path.join(r'result\urban\classification', 'gt_map.png'),format='png',dpi=600)#bbox_inches='tight',pad_inches=0)
plt.close()
print('decode map get finished')