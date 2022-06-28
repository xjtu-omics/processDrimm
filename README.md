# 获得IAGS input (windows)

IAGS所使用的block文件可通过DRIMM-synteny生成。

[English version](./README_en.md)

## 1.生成DRIMM-synteny的input文件 (processOrthofinder.py)

```python
dir = './example'
sp = ['Brachy','Maize','Rice','Sorghum','Telongatum']
gff_list = ['Brachy.gff','Maize.gff','Rice.gff','Sorghum.gff','Telongatum.gff']
sp_ratio = [2,4,2,2,2]
```

+ dir：所有物种的gff文件存放路径(gff文件与[MCScanX](https://github.com/wyp1125/MCScanx) 输入相同)
+ sp：物种列表
+ Orthogroups.tsv：[OrthoFinder](https://github.com/davidemms/OrthoFinder) 的输出文件
+ sp_ratio：目标物种拷贝数，例如没有加倍事件为1，一次WGD事件为2

processOrthofinder生成三种文件，分别为：  
+ .all.sequence：完整物种同源基因序列。每一行表示一条染色体，基因以同源基因ID表示
+ .all.sequence.genename：与.all.sequence格式一致，但是基因以基因名表示
+ .sequence：过滤后物种同源基因序列。根据物种目标拷贝数对.all.sequence中超过目标拷贝数的同源基因进行过滤，以同源基因ID表示

随后将所有的.sequence文件进行合并，得到drimm.sequence文件作为DRIMM-synteny的输入。

## 2.运行DRIMM-synteny (drimm.exe)

+ 在./drimm路径下使用`dotnet publish -c Release -r win-x64`进行编译得到drimm.exe(生成于./drimm/drimm/bin/Release/netcoreapp3.1/win-x64)

运行drimm.exe

<div align="center">
<img src="https://s1.328888.xyz/2022/04/23/2R9HB.png">
</div>

DRIMM-synteny的输入：
+ drimm.sequence文件路径
+ 输出目录
+ cycleLengthThreshold控制共线块的连续性(默认值为20)
+ dustThreshold控制最大同源基因多样性，会将同源基因个数中超过这个上限的同源基因组进行过滤。对于上面的例子，各个物种的目标拷贝数为2,4,2,2,2，所以dustThreshold设置为13(2+4+2+2+2+1)

DRIMM-synteny的输出：
+ synteny.txt

```
0:1 2147483647 -1 
1:1 8997 
2:1 -2147482956 -2147482957 -2147482958 -2147482959 -2147482960 -2147482961 -2147482962 -2147482963 -2147482964 -2147482965 -2147482966 -2147482967 -2147482968 -2147482969 -2147482970 -2147482971 -2147482972 -2147482973 -2147482974 -2147482975 -2147482976 
3:1 13651 10527 2147481608 
4:1 10498 6507 
5:1 7035 2147482195 3426 -2147482913 -2147482912 -2147482911 -2147482910 -2147482909 -2147482908 -2147482907 -2147482906 -2147482905 -2147482904 -2147482903 -2147482902 -2147482901 -2147482900 -2147482899 -2147482898 -2147482897 -2147482896 -2147482895 -2147482894 -2147482893 
6:1 10464 
7:1 2147480311 
8:1 8926 8925 
9:1 6472 2147482102 992 
10:1 2147481278 13638 2147479645 128 
11:1 8808 13137 
12:1 2147481394 -2147483438 -2147483437 -2147483436 -2147483435 -2147483434 -2147483433 -2147483432 -2147483431 -2147483430 -2147483429 -2147483428 -2147483427 -2147483426 -2147483425 -2147483424 -2147483423 -2147483422 -2147483421 -2147483420 -2147483419 -2147483418 2147483493 10949 
13:1 2147478991 

# 第一列为block的ID
# 第二列(冒号后)为block在全部物种中出现次数
# 第三列及以后为block内部同源基因的ID
```
+ blocks.txt

```
-261 -434 -774 -1204 -1358 -1227 -1393 479 -1227 -1393 -1031 -1159 676 -971 -561 921 911 782 -918 1317 -827 -518 -897 1324 -1381 1311 -1377 -1310 -851 300 -860 -1287 -1438 -1430 1421 -843 -430 381 668 421 321 539 415 -472 -646 -475 -1001 -1303 -1072 -796 -1075 -1081 -1271 1085 -1091 1272 1051 -1282 -1379 1291 -1150 629 -1173 -934 -1335 -958 -923 1220 -1066 -1062 -768 1206 -607 -1054 -1230 -1048 755 -1044 -729 -1007 -658 -998 -815 -1035 -605 -638 -1023 -1012 1056 -713 -1196 -1176 686 -1100 -662 -981 -368 -1200 -285 -1125 -777 -1170 -641 -865 -1223 -848 -687 -951 -681 -950 -943 -714 -588 -698 -966 -647 -957 -810 -576 -899 880 -731 885 -580 -920 -914 -913 -912 -570 764 -1437 -783 916 -776 -1437 -783 916 590 1258 -878 -832 -792 -230 -348 -890 1225 -926 -1332 963 1328 936 -1323 -975 574 423 -245 -45 226 44 
-716 937 -571 -942 -1312 -944 -602 -947 -675 -949 -876 1373 905 -1337 -856 -665 622 603 -855 337 -873 -594 -846 -864 -733 -730 -1380 -1369 -869 -601 -215 733 -850 -980 769 -1256 -779 575 -1380 -1369 -450 -840 -1367 -1126 -797 -1128 1130 459 -1405 -1117 1187 706 -1195 717 -584 -398 1253 -586 1342 -1156 -591 696 -1408 1357 -1314 1039 1359 -1268 316 320 -417 -313 473 -653 403 -645 1422 -990 1030 -1302 -1006 -1040 999 1289 1397 -1288 -1078 -1286 395 965 454 -1089 -745 -1090 -742 -1093 790 -1095 -1067 -1045 -572 784 941 -781 -1053 -789 -1058 -462 -1063 1232 1064 1350 -1402 -1285 -1284 -1065 -762 -1061 -1060 -766 -1059 -1057 -771 -772 1445 -1415 -1392 -1055 -1280 1392 1415 788 -788 -1445 -1047 -1277 -436 597 -1276 -734 -739 -1092 -482 -480 1088 -748 -750 -1087 -1273 -1067 1350 -1402 -255 925 -16 1095 -790 1093 742 1090 745 1089 -454 -965 -1286 380 1078 1288 -1397 -1289 -999 1040 1006 1302 1030 990 -1422 645 199 1388 290 653 17 
-1024 -1304 1004 1203 -239 1041 -997 496 -812 -995 -994 -993 -821 1298 -1297 -991 988 -1296 -986 -499 -985 461 749 1008 -723 -1009 -1011 829 583 1187 706 -1195 717 -584 -398 1253 -586 1342 204 474 1422 696 -1408 1357 -1314 1039 1359 -1268 -643 -655 1038 1290 1037 -650 582 1099 -630 1022 996 806 -540 820 -1019 -625 -1388 -1017 1344 1016 1216 956 1321 1013 -617 1020 287 1076 -286 1129 1174 1189 -1245 -1363 -1166 -1164 801 -1161 -1397 -1157 -722 -1154 1348 1248 -1175 -1178 -715 1172 -716 937 571 -942 -1312 -944 675 947 -949 -876 200 873 440 855 622 603 665 856 1337 -905 -1373 431 -594 -846 840 733 -850 -864 -601 -730 -1380 -1369 -869 455 -575 779 1256 -769 980 -1367 -1126 -797 -1128 537 1130 -597 389 -1276 -734 -739 -1092 480 482 1088 -748 -750 -1087 -1273 -925 -1183 874 1235 -688 -674 -1149 -680 394 274 
-1250 1116 1327 -1115 -791 -1113 -1266 -1112 -679 1111 1265 -22 1217 -305 1231 1106 -1105 -1262 1005 719 1073 847 1050 1046 1251 -1101 -1375 -1122 -1307 -1146 1259 -697 -1145 -657 -644 1182 628 627 -1257 -1141 -1139 -564 -1137 807 -566 1136 1267 -1132 1293 -535 798 -419 844 -830 -754 -1103 -1107 -743 -741 1110 -1120 -1179 -1441 -1152 868 866 894 1281 297 -21 
-21 -1357 1408 770 1165 -894 -866 1441 868 1152 1441 1179 -1120 -1110 741 743 1107 505 -409 844 -830 -754 -1103 505 -409 577 311 533 577 -535 798 -1293 1132 -1267 -1136 807 -566 1137 -608 -1139 -564 -608 1139 1141 1257 -627 -628 -1182 644 -697 -1145 -657 -1259 1146 1307 1122 1375 1101 -1251 -1046 -1050 467 -1020 617 -1013 -1321 -956 -1216 -1016 -1344 1017 1388 625 1019 -820 540 -806 -996 -203 -1022 630 -1099 1099 -1037 -1290 -1038 655 643 -1172 715 1178 1175 -1248 -1348 1154 -722 1157 1189 -631 -611 -862 1032 -689 606 -700 -1438 -1430 1421 490 843 -1421 1430 1438 1287 860 445 851 1310 1377 -1311 1381 -1324 897 238 921 911 782 -918 1317 971 -676 1159 1031 1393 1227 1358 1204 -1393 479 1358 1204 774 434 261 
-423 -574 975 1323 -936 -1328 -963 1332 926 -1225 890 -842 -590 418 -740 1437 792 832 878 -1258 -590 -916 783 1437 -764 570 912 -914 -913 920 880 -731 885 193 899 576 -588 698 -413 -966 -647 -957 -810 714 -848 950 681 951 -943 579 1223 865 -641 1170 777 1125 593 -1200 981 662 1100 -686 1176 1196 713 -1056 1012 1023 -301 -605 -638 729 -658 -998 815 -1035 -1007 1044 -755 1048 1230 1054 607 -1206 768 1062 1066 -1220 923 958 1335 934 1173 -629 1150 -1291 1379 1282 -1051 -1272 1091 -1085 1271 1081 1075 796 1072 1303 1001 -415 -539 -321 -421 -668 -381 430 
-933 -932 721 -704 -702 -940 -1319 -954 -695 -694 955 -1333 -976 -616 -974 -623 849 1162 1322 1079 973 -528 -875 1180 -621 -1353 -618 -969 893 -527 -1279 928 838 364 -526 854 523 -420 -189 -442 545 1398 -959 -362 -1398 1080 1077 1002 -339 -929 1094 898 1238 -1427 -1331 726 -895 1097 -517 -886 -1345 1383 642 889 881 879 1131 -673 -887 -209 900 809 -902 1123 922 -1346 -919 -917 -1340 -220 -1297 -991 988 -1296 -986 -1011 985 1009 723 -812 -361 -583 -829 461 749 1008 997 -234 994 995 -821 1298 867 549 463 786 232 -699 -904 243 -1041 492 -632 -1203 -1004 1304 1024 800 516 898 -426 282 1240 
93 -186 -393 -562 -1250 1116 1327 -1115 -791 -1113 -1266 -1112 -679 1111 1265 1217 560 1231 1106 -1105 -1262 1005 719 1073 847 1076 514 1129 318 1174 1161 -801 1164 1166 1363 1245 -841 -559 -336 -1363 -198 -811 932 933 721 -704 695 954 1319 -955 694 -702 -940 -1333 -976 -616 -974 -623 849 1162 1322 1079 973 969 -149 -1180 875 893 838 -1279 -527 -1279 928 526 456 319 -382 -786 -463 -549 -867 1340 917 919 1346 -922 1123 900 809 -902 511 1131 -673 -887 509 -879 -881 -889 -642 -1383 1345 -263 -1383 1345 886 347 258 345 1097 -345 895 -726 1331 1427 -1238 633 1034 443 344 727 -968 1329 -972 1010 

# 物种的共线块信息，其中正负号表示方向
```

## 3.生成IAGS的input (processDrimm.py)

```python
block_file = './example/drimm/blocks.txt'
drimmSyntenyFile = './example/drimm/synteny.txt'
outdir = './example'
chr_number = [5,10,12,10,7]
sp_list = ['Brachy','Maize','Rice','Sorghum','Telongatum']
target_rate = '2:4:2:2:2'
```
+ block_file：上一步生成的blocks.txt路径
+ drimmSyntenyFile：上一步生成的synteny.txt路径
+ outdir：结果文件输出路径(应与processOrthofinder.py中的输出路径一致)
+ chr_number：各个物种的染色体数
+ sp_list：物种列表
+ target_rate：各目标拷贝数比例

processDrimm首先将blocks文件进行拆分，得到各个物种的.block文件，并将.block文件中的共线块序列还原为同源基因ID序列，由于processOrthofinder步骤中过滤了不符合拷贝数的基因，故使用LCS算法将其和物种实际同源基因ID序列进行匹配，恢复各个共线块在原始物种中的同源基因ID序列以及genename序列，可供下游生物学分析。最后，将各物种.block文件中对不满足期望拷贝数的共线块过滤，得到各个物种的.final.block，供IAGS输入。  


+ LCS示意图
  ```
  1720:1 24 1256
  1640:1 33 666 1230 720 140
  666:1 557 886
  
  # blocks.txt中共线块ID和同源基因的对应关系
  ```

  ![LCS matching](https://s1.328888.xyz/2022/06/16/pxYEZ.png)

  <div align=center>图1 共线块完整同源基因序列的构建。</div>

  ```
  1720:1 24 288 1256
  1640:1 366 666 1230 720 2348 140
  666:1 557 567 886
  
  # LCS处理后，共线块在物种A实际的同源基因ID
  ----------------------------------------------------
  1720:1 a1 a2 a3
  1640:1 a4 a5 a6 a7 a8 a9
  666:1 a10 a11 a12
  
  # LCS处理后，共线块在物种A实际的同源基因名
  ```
  
  ```
  1720:1 24 357 1256
  1640:1 666 1230 543 720 1440 140
  666:1 557 886
  
  # LCS处理后，共线块在物种B实际的同源基因ID
  ----------------------------------------------------
  1720:1 b1 b2 b3
  1640:1 b4 b5 b6 b7 b8 b9
  666:1 b11 b12
  
  # LCS处理后，共线块在物种B实际的同源基因名
  ```

最终输出drimmBlocks和finalBlocks两个文件夹，两个文件夹中文件格式是相同的。不同点在于，drimmBlocks中的.block文件为未经过滤的block序列。另外两个.synteny和.synteny.genename文件为这个.block文件对应的共线块信息。finalBlocks中的.final.block文件为过滤后符合物种期望拷贝数的block，而另外两种.synteny和.synteny.genename文件为其.final.block文件对应的共线块信息。其中.final.block用作IAGS的输入文件。


文件格式示例：

+ Brachy.block

```
s -670 -1019 -1360 -1025 -1369 -632 -1344 360 -571 -736 1351 -1337 -1349 -852 -1368 -1182 739 -1629 -732 -1590 -768 -721 1487 847 -771 252 -777 -817 -1592 -476 -498 -1239 -1220 -1378 -1388 1363 -1622 -752 -1649 1418 1580 -1062 -1635 -1481 -1358 620 -1073 -583 -435 315 -354 881 -592 -343 -975 -450 849 -1560 -988 -1561 -999 -613 -1390 -957 1384 1679 -993 -1192 -904 -1482 887 -1485 -1343 1309 -1281 1147 909 1304 1118 1218 -1604 -1376 -1288 1376 1604 1468 -1289 -1291 -1058 542 -1063 972 -1018 700 -471 692 -324 -687 1069 -387 -81 451 1021 1495 1037 291 -1088 1401 -1080 -1467 1285 -1479 1554 -1462 -1081 1366 1240 1562 -622 -1091 1272 913 -1662 1631 -1540 -1339 -931 718 -1341 -679 -1009 286 -1539 -1483 -1565 1341 -882 -1398 -1488 -1526 -1396 -1692 -1162 -889 645 -1512 -1448 1492 -614 -1319 -944 -1447 -651 -646 1559 -639 -1321 -1515 -1098 -950 -625 -953 -1505 -954 -1405 1406 -963 1408 843 -1564 -1409 1310 1348 -1535 -973 -1538 -1508 -1546 -968 -967 1544 1542 -964 -595 -960 -1541 -958 -603 -1383 -1543 1392 -606 -1359 -948 -947 -946 -945 -100 1094 892 -977 77 995 704 -1360 626 1093 -1004 -1369 -1003 -691 829 695 1351 -1337 -1349 -709 -1368 -1182 -91 -438 448 382 994 -682 1590 -989 1629 -986 -985 1000 -984 1487 -983 1289 -1468 -1604 -1376 -1288 -1256 -1304 -909 -1147 1281 -1309 1343 1485 -538 1482 915 518 -1679 -1384 600 1390 -1153 -902 1561 -901 1560 -900 766 1154 -896 -1316 -11 996 894 -516 714 -1333 -1312 -1216 730 -1335 908 -497 1466 -1529 1166 -941 1442 -939 -938 -1506 -1495 -936 1443 -1504 -1501 -1446 -930 -1502 -1243 -1510 -1241 401 -1524 1024 -1522 978 1521 916 2 
s -912 1164 -942 -1010 -1227 1090 -828 1450 1020 684 -1074 1194 -1416 -1519 -1518 1480 -535 -1206 -1255 -1231 -1417 -1395 -550 -1039 -1221 -1035 -1030 827 924 547 -1426 -543 -541 407 -540 -1425 1305 1424 -536 -492 -1138 -1423 1300 -1422 1299 1016 1174 -556 -1428 -410 -587 1705 -1012 -1138 -587 1705 -1012 -524 1038 1168 1633 1632 974 -1067 1551 -1056 -1630 -1052 -320 -1146 1523 876 -1628 -1045 -1099 -509 -1429 -1626 -1042 -1625 -1431 -1624 -1445 -1444 -1672 -1620 -744 -1619 -734 -729 -1641 -725 -1644 -1566 -724 -1643 -746 1557 747 -767 1652 -765 -1648 -722 1514 -1449 -1665 -1591 -1664 -1150 -1663 -1172 -1435 1308 -1647 1434 -1650 -1407 -1653 -1618 -842 -1180 -841 840 1414 -1402 1528 -1393 -545 -831 -890 -1424 -1305 1425 -830 878 228 -1257 -534 530 1214 1395 1417 -866 -1114 -1480 1518 1519 1416 -531 -864 -1423 1300 -1422 1299 -584 -857 -1428 1022 1585 914 -1426 1617 970 1028 1071 1077 -1582 1096 578 570 336 428 -794 -1112 -1705 -787 1017 -562 -1040 -561 1633 1041 1632 1084 -780 1293 1066 1551 1054 -1630 779 -778 -1143 -1146 1523 -783 -1628 838 -1429 -1626 -797 -1625 -1431 -1624 -1445 -1444 -1672 -1620 -798 792 -1619 -821 1117 -820 -1641 -818 -1644 -1566 776 -1643 -785 1557 891 737 1652 764 -1648 -816 1514 -1449 -1665 -1591 -1664 -1150 -1663 -808 -1435 1308 -1647 1434 -1650 -1407 -1653 -1618 965 1013 1393 -1528 1402 -1414 -802 -1061 1582 1078 -1382 1617 -801 
s 181 -803 920 -671 520 -811 880 -660 -473 1339 1540 -1631 1662 -1179 -1272 935 -1562 -1032 1462 -1554 1479 -1285 1467 -1401 623 1043 -1322 -1539 -1483 -1565 -774 -1398 -1488 -1526 -1396 -1692 -1203 -1512 -1448 1492 -654 -1319 481 -1447 1321 567 1170 -1559 1155 -1515 617 420 -1160 -580 44 329 1228 489 1244 1205 -599 1585 1085 851 966 -1097 758 -1235 -806 -1608 -860 -861 1705 1083 1007 -862 246 616 -1489 1611 987 529 -436 -865 -480 -636 893 -805 -1119 -371 1450 -1264 -1667 272 -1282 1316 -789 -1207 1210 -836 1246 1187 1667 491 -1333 -1312 314 657 1466 -1529 -664 -845 1442 824 1524 1026 1241 519 -675 1510 1243 1502 1443 -1504 -1501 -1446 1059 1506 -823 -1522 917 1521 756 1001 699 710 -1123 -693 -621 -1453 -602 -1454 1345 1456 -711 -1457 1459 680 751 -404 -1460 -751 -643 1411 -1301 -1248 -1615 1575 1463 -1609 -1205 -1244 -604 1413 1520 1095 -1577 -1579 -1581 -733 1594 -1302 1572 -726 1589 -1375 -1656 -731 1470 -1471 1472 -411 -668 1031 863 -489 1121 -1269 -1185 -1181 1433 -1101 -1639 1377 551 1332 1412 -673 -1347 527 -1356 1614 813 -1651 -738 -665 1440 -727 -1588 839 1645 -762 -1578 898 1587 1473 -1573 1593 1612 -598 1576 -757 -648 -1602 -749 1598 -723 716 -1600 -719 1051 -1354 -720 1584 952 -1505 759 -609 -1405 772 -1541 815 -1383 -1543 1392 -846 -1359 859 1027 -1408 -844 -1406 -835 1409 1564 -832 728 -1542 -1544 -872 969 1546 1508 1538 -870 1535 -1348 -1310 959 761 479 -130 
s 804 -1592 1189 393 -1378 -1388 1363 -1622 754 -775 -1649 1418 1580 -991 -1635 -1481 -1358 1011 1605 -800 781 1627 -1278 -1159 -807 -812 -1141 475 -514 -350 1371 357 669 -263 -416 853 -1667 437 -1207 1210 -466 389 546 -1072 -1605 784 -826 1627 418 -434 -1240 -1366 1278 -1627 826 662 -769 703 507 -1173 -383 -485 -1124 624 619 -513 998 610 446 -356 392 444 1065 -502 -511 -1044 1251 1055 -1047 -1060 -422 1313 -1374 -1634 1372 1509 -1427 1295 -1595 1661 1511 1438 1637 1029 -1700 -1693 1513 1365 -449 686 -328 -1187 -1246 1103 -1218 -1118 -921 1123 715 -1453 -925 -1454 1345 1456 -711 -1457 1459 -577 1460 -572 956 1411 -1301 -1248 -1615 753 1575 1463 -1609 850 1413 1520 -933 -1577 -1579 -1581 -940 1594 -1302 1572 -911 1589 -1375 -1656 -884 1470 -1471 1472 1107 -1280 -1163 1075 -895 -899 1008 -368 427 303 -394 -442 -1371 -1060 -400 -649 1313 -1374 -1634 1372 1509 -1427 1295 -1595 1661 1294 1693 1700 1260 -1637 -1438 -1511 1184 -1329 1247 1513 1365 -1191 219 
s -219 589 1049 1129 992 982 745 310 -628 -949 -403 681 1228 448 382 -1181 1433 742 -1639 1377 1048 1332 1412 1057 -1347 712 -1356 926 -1614 -888 -1651 928 1440 717 -1588 822 1645 -885 -1578 -883 740 1587 1473 -1573 1593 1612 -1158 1576 -934 -1602 1079 1598 -927 559 1036 -1600 1064 683 -923 -814 -697 533 -1354 886 1584 1046 -1611 1489 -1023 -1076 -1235 1092 -1608 760 867 -413 591 -1086 -825 501 440 

# 物种的共线块信息，其中正负号表示方向，s表示线性染色体
```
+ Brachy.final.block
```
s -1590 -1592 -1622 -1649 -1635 -1561 1554 1562 -1662 1631 -1540 -1565 -1692 -1512 -1505 -1564 -1535 -1546 1544 1542 -1541 1590 1561 -1506 -1502 -1510 -1524 -1522 1521 
s -1519 -1518 1632 1551 -1630 1523 -1628 -1626 -1625 -1624 -1619 -1643 1557 1652 -1648 1514 -1665 -1664 -1663 -1653 -1618 1528 1518 1519 -1582 1632 1551 -1630 1523 -1628 -1626 -1625 -1624 -1619 -1643 1557 1652 -1648 1514 -1665 -1664 -1663 -1653 -1618 -1528 1582 
s 1540 -1631 1662 -1562 -1554 -1565 -1692 -1512 1524 1510 1502 1506 -1522 1521 1520 -1577 -1579 1594 1589 -1656 -1651 -1588 1645 -1578 1587 -1573 1593 1612 -1602 1598 -1600 -1505 -1541 1564 -1542 -1544 1546 1535 
s -1592 -1622 -1649 -1635 -1634 1509 1511 -1700 1520 -1577 -1579 1594 1589 -1656 -1634 1509 1700 -1511 
s -1651 -1588 1645 -1578 1587 -1573 1593 1612 -1602 1598 -1600 

# 物种的共线块信息，其中正负号表示方向，s表示线性染色体
# .final.block中记录的是过滤后符合物种期望拷贝数的block，可用作IAGS输入，格式与drimmBlocks文件夹中的.block相同
```

+ Brachy.synteny
```
670:1:chr_1 8450 8449 8448 18178 8447 19377 4760 16508 8446 19376 16507 8445 
1019:1:chr_1 8459 18180 8458 8457 8456 8455 4764 19378 2732 4763 8454 18179 8453 8452 3513 1492 8451 4762 153 4761 
1360:1:chr_1 2057 2056 16509 18181 8460 3514 
1025:1:chr_1 8505 18185 19379 18185 8504 8503 3523 8502 8501 4786 8500 4785 8499 8498 8497 903 8496 8495 8494 3522 4784 2058 4783 4782 8493 4781 16516 4780 8492 4779 8491 4778 651 651 8490 4777 8489 1494 4776 8488 8487 8486 3521 8485 3520 651 651 1493 4775 3519 8484 3518 16515 4774 4773 2735 486 2734 2734 2734 8483 16514 4772 8482 4771 16513 8481 18184 8480 18184 8479 49 49 49 49 18183 8478 8477 16512 8476 4770 8475 8474 49 187 8473 8472 8471 2733 3517 902 3516 8470 18182 2733 4769 8469 4768 8468 4767 8467 8466 16511 4766 8465 8464 8463 4765 8462 16510 3515 8461 
1369:1:chr_1 2737 551 25 19381 8518 16520 1497 8517 8516 8515 8514 16519 2060 8513 8512 4788 8511 19380 3524 8510 1496 2736 8509 769 8508 8507 1495 8506 4787 
1369:2:chr_1 2737 215 215 16740 5524 5395 2 584 16739 2223 2223 4791 10082 551 19468 10081 3631 10080 16738 10079 10078 1495 

# drimmBlocks文件夹中，Brachy.block对应的block信息(同源基因ID形式)
# 第一列为block的ID
# 第二列(第一个冒号后)为block在物种中出现的序数，计数顺序为从左到右，从上到下
# 第三列(第二个冒号后)为block所在的染色体号
# 第四列(第一个空格后)及以后为LCS匹配后，block在实际物种的同源基因的ID
```
+ Brachy.synteny.genename
```
670:1:chr_1 BrachyLOC100841065 BrachyLOC100826686 BrachyLOC106865846 BrachyLOC104581777 BrachyLOC100838523 BrachyLOC100837927 BrachyLOC100821843 BrachyLOC100844753 BrachyLOC100837316 BrachyLOC100837010 BrachyLOC100836696 BrachyLOC100825490 
1019:1:chr_1 BrachyLOC100843403 BrachyLOC100843101 BrachyLOC100825346 BrachyLOC100833933 BrachyLOC100840563 BrachyLOC100844407 BrachyLOC100842797 BrachyLOC100832900 BrachyLOC100837517 BrachyLOC100842179 BrachyLOC104582101 BrachyLOC100841885 BrachyLOC100829233 BrachyLOC100832494 BrachyLOC100841577 BrachyLOC100822991 BrachyLOC100841273 BrachyLOC100823300 BrachyLOC100844504 BrachyLOC100840970 
1360:1:chr_1 BrachyLOC100844822 BrachyLOC104582482 BrachyLOC100825655 BrachyLOC100843713 BrachyLOC100841683 BrachyLOC100833436 
1025:1:chr_1 BrachyLOC100838735 BrachyLOC106865425 BrachyLOC112268512 BrachyLOC100838129 BrachyLOC100837211 BrachyLOC100839980 BrachyLOC100836616 BrachyLOC100833662 BrachyLOC100834883 BrachyLOC100838363 BrachyLOC100828463 BrachyLOC100826418 BrachyLOC100829986 BrachyLOC100824882 BrachyLOC100836288 BrachyLOC106866865 BrachyLOC100844540 BrachyLOC100840394 BrachyLOC100837455 BrachyLOC100834746 BrachyLOC100831008 BrachyLOC100834442 BrachyLOC100834142 BrachyLOC100833833 BrachyLOC100833526 BrachyLOC100833215 BrachyLOC104583610 BrachyLOC100832904 BrachyLOC100824159 BrachyLOC100822097 BrachyLOC100846376 BrachyLOC104583539 BrachyLOC104581932 BrachyLOC104581931 BrachyLOC100843234 BrachyLOC100841811 BrachyLOC100837453 BrachyLOC100834382 BrachyLOC100832282 BrachyLOC100835607 BrachyLOC100828570 BrachyLOC100831924 BrachyLOC100826200 BrachyLOC100823545 BrachyLOC100831992 BrachyLOC104581930 BrachyLOC104581929 BrachyLOC100843844 BrachyLOC100846578 BrachyLOC100831385 BrachyLOC100831082 BrachyLOC100838562 BrachyLOC100836116 BrachyLOC100829675 BrachyLOC100830155 BrachyLOC100827039 BrachyLOC100829849 BrachyLOC112270139 BrachyLOC100828940 BrachyLOC104581927 BrachyLOC100839483 BrachyLOC100828638 BrachyLOC100828333 BrachyLOC100828030 BrachyLOC100832301 BrachyLOC100827719 BrachyLOC112268584 BrachyLOC100827109 BrachyLOC100823025 BrachyLOC100838057 BrachyLOC100825068 BrachyLOC100825552 BrachyLOC112268527 BrachyLOC100825248 BrachyLOC100824945 BrachyLOC100824336 BrachyLOC100823118 BrachyLOC100833957 BrachyLOC100842093 BrachyLOC100822796 BrachyLOC100834672 BrachyLOC100843616 BrachyLOC100824036 BrachyLOC100838051 BrachyLOC104582792 BrachyLOC100836408 BrachyLOC100822170 BrachyLOC100839468 BrachyLOC100824755 BrachyLOC100821862 BrachyLOC100821551 BrachyLOC104581926 BrachyLOC100845650 BrachyLOC100846757 BrachyLOC100846455 BrachyLOC100846148 BrachyLOC100829249 BrachyLOC100838851 BrachyLOC100846048 BrachyLOC100826708 BrachyLOC100836597 BrachyLOC100844129 BrachyLOC100845846 BrachyLOC100834363 BrachyLOC100844620 BrachyLOC100830366 BrachyLOC100825155 BrachyLOC100835781 BrachyLOC100826703 BrachyLOC100844316 BrachyLOC100844018 BrachyLOC100837430 
1369:1:chr_1 BrachyLOC100846038 BrachyLOC100833360 BrachyLOC100845732 BrachyLOC100845434 BrachyLOC100826420 BrachyLOC100822411 BrachyLOC112270045 BrachyLOC100841612 BrachyLOC100844212 BrachyLOC100832620 BrachyLOC100837657 BrachyLOC100834472 BrachyLOC100843911 BrachyLOC100828669 BrachyLOC100826204 BrachyLOC100823739 BrachyLOC100821579 BrachyLOC112270150 BrachyLOC100843602 BrachyLOC100843301 BrachyLOC100842996 BrachyLOC100836738 BrachyLOC100842692 BrachyLOC100842386 BrachyLOC100842077 BrachyLOC100824160 BrachyLOC104584204 BrachyLOC100840864 BrachyLOC100840704 
1369:2:chr_1 BrachyLOC100841064 BrachyLOC100840760 BrachyLOC100840150 BrachyLOC100839851 BrachyLOC100839548 BrachyLOC112270598 BrachyLOC100829059 BrachyLOC106865684 BrachyLOC100828757 BrachyLOC100839242 BrachyLOC100838935 BrachyLOC100838625 BrachyLOC100828154 BrachyLOC100838115 BrachyLOC100827544 BrachyLOC100827237 BrachyLOC100837811 BrachyLOC100837508 BrachyLOC100826619 BrachyLOC100837201 BrachyLOC104582271 BrachyLOC100836581 

# drimmBlocks文件夹中，Brachy.block对应的block信息(基因名形式)
# 第一列为block的ID
# 第二列(第一个冒号后)为block在物种中出现的序数，计数顺序为从左到右，从上到下
# 第三列(第二个冒号后)为block所在的染色体号
# 第四列(第一个空格后)及以后为LCS匹配后，block在实际物种中的基因名
```


## 4.生成IAGS绘图所需的blockindex.genenumber文件

```python
sp = ['Brachy','Maize','Rice','Sorghum','Telongatum']
finalBlocks_path = './example/finalBlocks'
```
+ sp：物种列表
+ finalBlocks_path：上一步生成的finalBlocks文件夹路径

processGenenumber的输出为blockindex.genenumber文件，文件记录了共线块号与其在实际物种中最大基因数量的对应关系。该文件可以用作IAGS画图功能中plotChrsRearrangement的block_length_file参数。

+ blockindex.genenumber

```
blockID	blockLength
1590	46
1592	50
1622	16
1649	114
1635	23
1561	36
1554	175

# 第一列为block号，第二列为block长度，其中block长度为LCS匹配后，block在实际物种中最大基因数量
```

## DRIMM-synteny

[DRIMM-Synteny: decomposing genomes into evolutionary conserved segments ](https://academic.oup.com/bioinformatics/article/26/20/2509/193644?login=false)


[DRIMM - Duplications and Rearrangements In Multiple Mammals](http://bix.ucsd.edu/projects/drimm/)
