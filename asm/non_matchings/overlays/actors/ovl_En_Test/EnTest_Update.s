glabel EnTest_Update
/* 04044 80863694 27BDFFB0 */  addiu   $sp, $sp, 0xFFB0           ## $sp = FFFFFFB0
/* 04048 80863698 AFBF002C */  sw      $ra, 0x002C($sp)           
/* 0404C 8086369C AFB00028 */  sw      $s0, 0x0028($sp)           
/* 04050 808636A0 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 04054 808636A4 0C218D3E */  jal     func_808634F8              
/* 04058 808636A8 AFA50054 */  sw      $a1, 0x0054($sp)           
/* 0405C 808636AC 920E00B1 */  lbu     $t6, 0x00B1($s0)           ## 000000B1
/* 04060 808636B0 24010006 */  addiu   $at, $zero, 0x0006         ## $at = 00000006
/* 04064 808636B4 51C100A7 */  beql    $t6, $at, .L80863954       
/* 04068 808636B8 26050810 */  addiu   $a1, $s0, 0x0810           ## $a1 = 00000810
/* 0406C 808636BC 0C00B638 */  jal     Actor_MoveForward
              
/* 04070 808636C0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 04074 808636C4 3C0141F0 */  lui     $at, 0x41F0                ## $at = 41F00000
/* 04078 808636C8 44810000 */  mtc1    $at, $f0                   ## $f0 = 30.00
/* 0407C 808636CC 240F001D */  addiu   $t7, $zero, 0x001D         ## $t7 = 0000001D
/* 04080 808636D0 AFAF0014 */  sw      $t7, 0x0014($sp)           
/* 04084 808636D4 44070000 */  mfc1    $a3, $f0                   
/* 04088 808636D8 8FA40054 */  lw      $a0, 0x0054($sp)           
/* 0408C 808636DC 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 04090 808636E0 3C064296 */  lui     $a2, 0x4296                ## $a2 = 42960000
/* 04094 808636E4 0C00B92D */  jal     func_8002E4B4              
/* 04098 808636E8 E7A00010 */  swc1    $f0, 0x0010($sp)           
/* 0409C 808636EC 8618001C */  lh      $t8, 0x001C($s0)           ## 0000001C
/* 040A0 808636F0 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 040A4 808636F4 5701001F */  bnel    $t8, $at, .L80863774       
/* 040A8 808636F8 96090088 */  lhu     $t1, 0x0088($s0)           ## 00000088
/* 040AC 808636FC C600000C */  lwc1    $f0, 0x000C($s0)           ## 0000000C
/* 040B0 80863700 C6040028 */  lwc1    $f4, 0x0028($s0)           ## 00000028
/* 040B4 80863704 4600203E */  c.le.s  $f4, $f0                   
/* 040B8 80863708 00000000 */  nop
/* 040BC 8086370C 45020005 */  bc1fl   .L80863724                 
/* 040C0 80863710 C6080080 */  lwc1    $f8, 0x0080($s0)           ## 00000080
/* 040C4 80863714 44803000 */  mtc1    $zero, $f6                 ## $f6 = 0.00
/* 040C8 80863718 E6000028 */  swc1    $f0, 0x0028($s0)           ## 00000028
/* 040CC 8086371C E6060060 */  swc1    $f6, 0x0060($s0)           ## 00000060
/* 040D0 80863720 C6080080 */  lwc1    $f8, 0x0080($s0)           ## 00000080
.L80863724:
/* 040D4 80863724 4600403E */  c.le.s  $f8, $f0                   
/* 040D8 80863728 00000000 */  nop
/* 040DC 8086372C 45020003 */  bc1fl   .L8086373C                 
/* 040E0 80863730 8E1907CC */  lw      $t9, 0x07CC($s0)           ## 000007CC
/* 040E4 80863734 E6000080 */  swc1    $f0, 0x0080($s0)           ## 00000080
.L80863738:
/* 040E8 80863738 8E1907CC */  lw      $t9, 0x07CC($s0)           ## 000007CC
.L8086373C:
/* 040EC 8086373C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 040F0 80863740 8FA50054 */  lw      $a1, 0x0054($sp)           
/* 040F4 80863744 0320F809 */  jalr    $ra, $t9                   
/* 040F8 80863748 00000000 */  nop
/* 040FC 8086374C 920207DE */  lbu     $v0, 0x07DE($s0)           ## 000007DE
/* 04100 80863750 2C410005 */  sltiu   $at, $v0, 0x0005           
/* 04104 80863754 10200067 */  beq     $at, $zero, .L808638F4     
/* 04108 80863758 00024080 */  sll     $t0, $v0,  2               
/* 0410C 8086375C 3C018086 */  lui     $at, %hi(jtbl_808647E4)       ## $at = 80860000
/* 04110 80863760 00280821 */  addu    $at, $at, $t0              
/* 04114 80863764 8C2847E4 */  lw      $t0, %lo(jtbl_808647E4)($at)  
/* 04118 80863768 01000008 */  jr      $t0                        
/* 0411C 8086376C 00000000 */  nop
/* 04120 80863770 96090088 */  lhu     $t1, 0x0088($s0)           ## 00000088
.L80863774:
/* 04124 80863774 8FA40054 */  lw      $a0, 0x0054($sp)           
/* 04128 80863778 312A0002 */  andi    $t2, $t1, 0x0002           ## $t2 = 00000000
/* 0412C 8086377C 1140FFEE */  beq     $t2, $zero, .L80863738     
/* 04130 80863780 248407C0 */  addiu   $a0, $a0, 0x07C0           ## $a0 = 000007C0
/* 04134 80863784 8E050078 */  lw      $a1, 0x0078($s0)           ## 00000078
/* 04138 80863788 9206007D */  lbu     $a2, 0x007D($s0)           ## 0000007D
/* 0413C 8086378C 0C0107A9 */  jal     func_80041EA4              
/* 04140 80863790 AFA40038 */  sw      $a0, 0x0038($sp)           
/* 04144 80863794 24010005 */  addiu   $at, $zero, 0x0005         ## $at = 00000005
/* 04148 80863798 1041000A */  beq     $v0, $at, .L808637C4       
/* 0414C 8086379C 8FA40038 */  lw      $a0, 0x0038($sp)           
/* 04150 808637A0 2401000C */  addiu   $at, $zero, 0x000C         ## $at = 0000000C
/* 04154 808637A4 10410007 */  beq     $v0, $at, .L808637C4       
/* 04158 808637A8 00000000 */  nop
/* 0415C 808637AC 8E050078 */  lw      $a1, 0x0078($s0)           ## 00000078
/* 04160 808637B0 0C010753 */  jal     func_80041D4C              
/* 04164 808637B4 9206007D */  lbu     $a2, 0x007D($s0)           ## 0000007D
/* 04168 808637B8 24010009 */  addiu   $at, $zero, 0x0009         ## $at = 00000009
/* 0416C 808637BC 5441FFDF */  bnel    $v0, $at, .L8086373C       
/* 04170 808637C0 8E1907CC */  lw      $t9, 0x07CC($s0)           ## 000007CC
.L808637C4:
/* 04174 808637C4 0C00B55C */  jal     Actor_Kill
              
/* 04178 808637C8 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0417C 808637CC 100000B6 */  beq     $zero, $zero, .L80863AA8   
/* 04180 808637D0 8FBF002C */  lw      $ra, 0x002C($sp)           
glabel L808637D4
/* 04184 808637D4 3C040600 */  lui     $a0, 0x0600                ## $a0 = 06000000
/* 04188 808637D8 0C028800 */  jal     SkelAnime_GetFrameCount
              
/* 0418C 808637DC 24841C20 */  addiu   $a0, $a0, 0x1C20           ## $a0 = 06001C20
/* 04190 808637E0 44825000 */  mtc1    $v0, $f10                  ## $f10 = 0.00
/* 04194 808637E4 3C014000 */  lui     $at, 0x4000                ## $at = 40000000
/* 04198 808637E8 44810000 */  mtc1    $at, $f0                   ## $f0 = 2.00
/* 0419C 808637EC 46805420 */  cvt.s.w $f16, $f10                 
/* 041A0 808637F0 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 041A4 808637F4 240B0002 */  addiu   $t3, $zero, 0x0002         ## $t3 = 00000002
/* 041A8 808637F8 44060000 */  mfc1    $a2, $f0                   
/* 041AC 808637FC AFAB0014 */  sw      $t3, 0x0014($sp)           
/* 041B0 80863800 24A51C20 */  addiu   $a1, $a1, 0x1C20           ## $a1 = 06001C20
/* 041B4 80863804 E7B00010 */  swc1    $f16, 0x0010($sp)          
/* 041B8 80863808 260404A8 */  addiu   $a0, $s0, 0x04A8           ## $a0 = 000004A8
/* 041BC 8086380C 24070000 */  addiu   $a3, $zero, 0x0000         ## $a3 = 00000000
/* 041C0 80863810 0C029468 */  jal     SkelAnime_ChangeAnimation
              
/* 041C4 80863814 E7A00018 */  swc1    $f0, 0x0018($sp)           
/* 041C8 80863818 3C0C8086 */  lui     $t4, %hi(D_80864510)       ## $t4 = 80860000
/* 041CC 8086381C 258C4510 */  addiu   $t4, $t4, %lo(D_80864510)  ## $t4 = 80864510
/* 041D0 80863820 92050188 */  lbu     $a1, 0x0188($s0)           ## 00000188
/* 041D4 80863824 8E0601A8 */  lw      $a2, 0x01A8($s0)           ## 000001A8
/* 041D8 80863828 8E0704C8 */  lw      $a3, 0x04C8($s0)           ## 000004C8
/* 041DC 8086382C AFAC0010 */  sw      $t4, 0x0010($sp)           
/* 041E0 80863830 0C028D52 */  jal     func_800A3548              
/* 041E4 80863834 8FA40054 */  lw      $a0, 0x0054($sp)           
/* 041E8 80863838 920D07DE */  lbu     $t5, 0x07DE($s0)           ## 000007DE
/* 041EC 8086383C 25AE0001 */  addiu   $t6, $t5, 0x0001           ## $t6 = 00000001
/* 041F0 80863840 1000002C */  beq     $zero, $zero, .L808638F4   
/* 041F4 80863844 A20E07DE */  sb      $t6, 0x07DE($s0)           ## 000007DE
glabel L80863848
/* 041F8 80863848 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 041FC 8086384C 260404A8 */  addiu   $a0, $s0, 0x04A8           ## $a0 = 000004A8
/* 04200 80863850 3C078086 */  lui     $a3, %hi(D_80864510)       ## $a3 = 80860000
/* 04204 80863854 24E74510 */  addiu   $a3, $a3, %lo(D_80864510)  ## $a3 = 80864510
/* 04208 80863858 26040188 */  addiu   $a0, $s0, 0x0188           ## $a0 = 00000188
/* 0420C 8086385C 8E0501A8 */  lw      $a1, 0x01A8($s0)           ## 000001A8
/* 04210 80863860 0C02950A */  jal     func_800A5428              
/* 04214 80863864 8E0604C8 */  lw      $a2, 0x04C8($s0)           ## 000004C8
/* 04218 80863868 10000023 */  beq     $zero, $zero, .L808638F8   
/* 0421C 8086386C 92180114 */  lbu     $t8, 0x0114($s0)           ## 00000114
glabel L80863870
/* 04220 80863870 3C014080 */  lui     $at, 0x4080                ## $at = 40800000
/* 04224 80863874 44819000 */  mtc1    $at, $f18                  ## $f18 = 4.00
/* 04228 80863878 244F0001 */  addiu   $t7, $v0, 0x0001           ## $t7 = 00000001
/* 0422C 8086387C A20F07DE */  sb      $t7, 0x07DE($s0)           ## 000007DE
/* 04230 80863880 E61204D0 */  swc1    $f18, 0x04D0($s0)          ## 000004D0
glabel L80863884
/* 04234 80863884 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 04238 80863888 44816000 */  mtc1    $at, $f12                  ## $f12 = 1.00
/* 0423C 8086388C C60004D0 */  lwc1    $f0, 0x04D0($s0)           ## 000004D0
/* 04240 80863890 44803000 */  mtc1    $zero, $f6                 ## $f6 = 0.00
/* 04244 80863894 460C0101 */  sub.s   $f4, $f0, $f12             
/* 04248 80863898 46000086 */  mov.s   $f2, $f0                   
/* 0424C 8086389C E60404D0 */  swc1    $f4, 0x04D0($s0)           ## 000004D0
/* 04250 808638A0 C60004D0 */  lwc1    $f0, 0x04D0($s0)           ## 000004D0
/* 04254 808638A4 4606003E */  c.le.s  $f0, $f6                   
/* 04258 808638A8 00000000 */  nop
/* 0425C 808638AC 45020004 */  bc1fl   .L808638C0                 
/* 04260 808638B0 46020203 */  div.s   $f8, $f0, $f2              
/* 04264 808638B4 A20007DE */  sb      $zero, 0x07DE($s0)         ## 000007DE
/* 04268 808638B8 C60004D0 */  lwc1    $f0, 0x04D0($s0)           ## 000004D0
/* 0426C 808638BC 46020203 */  div.s   $f8, $f0, $f2              
.L808638C0:
/* 04270 808638C0 8E0504C8 */  lw      $a1, 0x04C8($s0)           ## 000004C8
/* 04274 808638C4 92040188 */  lbu     $a0, 0x0188($s0)           ## 00000188
/* 04278 808638C8 8E0701A8 */  lw      $a3, 0x01A8($s0)           ## 000001A8
/* 0427C 808638CC 00A03025 */  or      $a2, $a1, $zero            ## $a2 = 00000000
/* 04280 808638D0 46086281 */  sub.s   $f10, $f12, $f8            
/* 04284 808638D4 0C028B9C */  jal     func_800A2E70              
/* 04288 808638D8 E7AA0010 */  swc1    $f10, 0x0010($sp)          
/* 0428C 808638DC 3C078086 */  lui     $a3, %hi(D_80864510)       ## $a3 = 80860000
/* 04290 808638E0 24E74510 */  addiu   $a3, $a3, %lo(D_80864510)  ## $a3 = 80864510
/* 04294 808638E4 26040188 */  addiu   $a0, $s0, 0x0188           ## $a0 = 00000188
/* 04298 808638E8 8E0501A8 */  lw      $a1, 0x01A8($s0)           ## 000001A8
/* 0429C 808638EC 0C02950A */  jal     func_800A5428              
/* 042A0 808638F0 8E0604C8 */  lw      $a2, 0x04C8($s0)           ## 000004C8
glabel L808638F4
.L808638F4:
/* 042A4 808638F4 92180114 */  lbu     $t8, 0x0114($s0)           ## 00000114
.L808638F8:
/* 042A8 808638F8 57000016 */  bnel    $t8, $zero, .L80863954     
/* 042AC 808638FC 26050810 */  addiu   $a1, $s0, 0x0810           ## $a1 = 00000810
/* 042B0 80863900 921900AF */  lbu     $t9, 0x00AF($s0)           ## 000000AF
/* 042B4 80863904 53200013 */  beql    $t9, $zero, .L80863954     
/* 042B8 80863908 26050810 */  addiu   $a1, $s0, 0x0810           ## $a1 = 00000810
/* 042BC 8086390C 920207C8 */  lbu     $v0, 0x07C8($s0)           ## 000007C8
/* 042C0 80863910 24010010 */  addiu   $at, $zero, 0x0010         ## $at = 00000010
/* 042C4 80863914 260407D2 */  addiu   $a0, $s0, 0x07D2           ## $a0 = 000007D2
/* 042C8 80863918 10410009 */  beq     $v0, $at, .L80863940       
/* 042CC 8086391C 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 042D0 80863920 24010017 */  addiu   $at, $zero, 0x0017         ## $at = 00000017
/* 042D4 80863924 50410007 */  beql    $v0, $at, .L80863944       
/* 042D8 80863928 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 042DC 8086392C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 042E0 80863930 0C218D18 */  jal     func_80863460              
/* 042E4 80863934 8FA50054 */  lw      $a1, 0x0054($sp)           
/* 042E8 80863938 10000006 */  beq     $zero, $zero, .L80863954   
/* 042EC 8086393C 26050810 */  addiu   $a1, $s0, 0x0810           ## $a1 = 00000810
.L80863940:
/* 042F0 80863940 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
.L80863944:
/* 042F4 80863944 240703E8 */  addiu   $a3, $zero, 0x03E8         ## $a3 = 000003E8
/* 042F8 80863948 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 042FC 8086394C AFA00010 */  sw      $zero, 0x0010($sp)         
/* 04300 80863950 26050810 */  addiu   $a1, $s0, 0x0810           ## $a1 = 00000810
.L80863954:
/* 04304 80863954 AFA50034 */  sw      $a1, 0x0034($sp)           
/* 04308 80863958 0C0189B7 */  jal     ActorCollider_Cylinder_Update
              
/* 0430C 8086395C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 04310 80863960 8E080028 */  lw      $t0, 0x0028($s0)           ## 00000028
/* 04314 80863964 3C014234 */  lui     $at, 0x4234                ## $at = 42340000
/* 04318 80863968 8E090024 */  lw      $t1, 0x0024($s0)           ## 00000024
/* 0431C 8086396C AE08003C */  sw      $t0, 0x003C($s0)           ## 0000003C
/* 04320 80863970 C610003C */  lwc1    $f16, 0x003C($s0)          ## 0000003C
/* 04324 80863974 44819000 */  mtc1    $at, $f18                  ## $f18 = 45.00
/* 04328 80863978 AE090038 */  sw      $t1, 0x0038($s0)           ## 00000038
/* 0432C 8086397C 8E09002C */  lw      $t1, 0x002C($s0)           ## 0000002C
/* 04330 80863980 46128100 */  add.s   $f4, $f16, $f18            
/* 04334 80863984 920A00AF */  lbu     $t2, 0x00AF($s0)           ## 000000AF
/* 04338 80863988 AE090040 */  sw      $t1, 0x0040($s0)           ## 00000040
/* 0433C 8086398C 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 04340 80863990 E604003C */  swc1    $f4, 0x003C($s0)           ## 0000003C
/* 04344 80863994 1D400003 */  bgtz    $t2, .L808639A4            
/* 04348 80863998 8FA40054 */  lw      $a0, 0x0054($sp)           
/* 0434C 8086399C 920B0114 */  lbu     $t3, 0x0114($s0)           ## 00000114
/* 04350 808639A0 1160001B */  beq     $t3, $zero, .L80863A10     
.L808639A4:
/* 04354 808639A4 34211E60 */  ori     $at, $at, 0x1E60           ## $at = 00011E60
/* 04358 808639A8 00812821 */  addu    $a1, $a0, $at              
/* 0435C 808639AC AFA50038 */  sw      $a1, 0x0038($sp)           
/* 04360 808639B0 0C017713 */  jal     Actor_CollisionCheck_SetOT
              ## CollisionCheck_setOT
/* 04364 808639B4 8FA60034 */  lw      $a2, 0x0034($sp)           
/* 04368 808639B8 920C07C8 */  lbu     $t4, 0x07C8($s0)           ## 000007C8
/* 0436C 808639BC 2981000A */  slti    $at, $t4, 0x000A           
/* 04370 808639C0 5420000D */  bnel    $at, $zero, .L808639F8     
/* 04374 808639C4 921807DE */  lbu     $t8, 0x07DE($s0)           ## 000007DE
/* 04378 808639C8 920D0114 */  lbu     $t5, 0x0114($s0)           ## 00000114
/* 0437C 808639CC 8FA40054 */  lw      $a0, 0x0054($sp)           
/* 04380 808639D0 8FA50038 */  lw      $a1, 0x0038($sp)           
/* 04384 808639D4 11A00005 */  beq     $t5, $zero, .L808639EC     
/* 04388 808639D8 00000000 */  nop
/* 0438C 808639DC 960E0112 */  lhu     $t6, 0x0112($s0)           ## 00000112
/* 04390 808639E0 31CF4000 */  andi    $t7, $t6, 0x4000           ## $t7 = 00000000
/* 04394 808639E4 55E00004 */  bnel    $t7, $zero, .L808639F8     
/* 04398 808639E8 921807DE */  lbu     $t8, 0x07DE($s0)           ## 000007DE
.L808639EC:
/* 0439C 808639EC 0C01767D */  jal     Actor_CollisionCheck_SetAC
              ## CollisionCheck_setAC
/* 043A0 808639F0 8FA60034 */  lw      $a2, 0x0034($sp)           
/* 043A4 808639F4 921807DE */  lbu     $t8, 0x07DE($s0)           ## 000007DE
.L808639F8:
/* 043A8 808639F8 8FA40054 */  lw      $a0, 0x0054($sp)           
/* 043AC 808639FC 8FA50038 */  lw      $a1, 0x0038($sp)           
/* 043B0 80863A00 53000004 */  beql    $t8, $zero, .L80863A14     
/* 043B4 80863A04 82190808 */  lb      $t9, 0x0808($s0)           ## 00000808
/* 043B8 80863A08 0C01767D */  jal     Actor_CollisionCheck_SetAC
              ## CollisionCheck_setAC
/* 043BC 80863A0C 260608DC */  addiu   $a2, $s0, 0x08DC           ## $a2 = 000008DC
.L80863A10:
/* 043C0 80863A10 82190808 */  lb      $t9, 0x0808($s0)           ## 00000808
.L80863A14:
/* 043C4 80863A14 5B200011 */  blezl   $t9, .L80863A5C            
/* 043C8 80863A18 860B001C */  lh      $t3, 0x001C($s0)           ## 0000001C
/* 043CC 80863A1C 9202086C */  lbu     $v0, 0x086C($s0)           ## 0000086C
/* 043D0 80863A20 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 043D4 80863A24 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 043D8 80863A28 30480004 */  andi    $t0, $v0, 0x0004           ## $t0 = 00000000
/* 043DC 80863A2C 15000008 */  bne     $t0, $zero, .L80863A50     
/* 043E0 80863A30 304AFFFB */  andi    $t2, $v0, 0xFFFB           ## $t2 = 00000000
/* 043E4 80863A34 8FA40054 */  lw      $a0, 0x0054($sp)           
/* 043E8 80863A38 34211E60 */  ori     $at, $at, 0x1E60           ## $at = 00011E60
/* 043EC 80863A3C 2606085C */  addiu   $a2, $s0, 0x085C           ## $a2 = 0000085C
/* 043F0 80863A40 0C0175E7 */  jal     Actor_CollisionCheck_SetAT
              ## CollisionCheck_setAT
/* 043F4 80863A44 00812821 */  addu    $a1, $a0, $at              
/* 043F8 80863A48 10000004 */  beq     $zero, $zero, .L80863A5C   
/* 043FC 80863A4C 860B001C */  lh      $t3, 0x001C($s0)           ## 0000001C
.L80863A50:
/* 04400 80863A50 0C218C90 */  jal     func_80863240              
/* 04404 80863A54 A20A086C */  sb      $t2, 0x086C($s0)           ## 0000086C
/* 04408 80863A58 860B001C */  lh      $t3, 0x001C($s0)           ## 0000001C
.L80863A5C:
/* 0440C 80863A5C 8FAC0054 */  lw      $t4, 0x0054($sp)           
/* 04410 80863A60 55600011 */  bnel    $t3, $zero, .L80863AA8     
/* 04414 80863A64 8FBF002C */  lw      $ra, 0x002C($sp)           
/* 04418 80863A68 918D1C27 */  lbu     $t5, 0x1C27($t4)           ## 00001C27
/* 0441C 80863A6C 3C188003 */  lui     $t8, 0x8003                ## $t8 = 80030000
/* 04420 80863A70 51A00008 */  beql    $t5, $zero, .L80863A94     
/* 04424 80863A74 8E190004 */  lw      $t9, 0x0004($s0)           ## 00000004
/* 04428 80863A78 8E0E0004 */  lw      $t6, 0x0004($s0)           ## 00000004
/* 0442C 80863A7C 2718B8C4 */  addiu   $t8, $t8, 0xB8C4           ## $t8 = 8002B8C4
/* 04430 80863A80 AE1800C0 */  sw      $t8, 0x00C0($s0)           ## 000000C0
/* 04434 80863A84 35CF0081 */  ori     $t7, $t6, 0x0081           ## $t7 = 00000081
/* 04438 80863A88 10000006 */  beq     $zero, $zero, .L80863AA4   
/* 0443C 80863A8C AE0F0004 */  sw      $t7, 0x0004($s0)           ## 00000004
/* 04440 80863A90 8E190004 */  lw      $t9, 0x0004($s0)           ## 00000004
.L80863A94:
/* 04444 80863A94 2401FF7E */  addiu   $at, $zero, 0xFF7E         ## $at = FFFFFF7E
/* 04448 80863A98 AE0000C0 */  sw      $zero, 0x00C0($s0)         ## 000000C0
/* 0444C 80863A9C 03214024 */  and     $t0, $t9, $at              
/* 04450 80863AA0 AE080004 */  sw      $t0, 0x0004($s0)           ## 00000004
.L80863AA4:
/* 04454 80863AA4 8FBF002C */  lw      $ra, 0x002C($sp)           
.L80863AA8:
/* 04458 80863AA8 8FB00028 */  lw      $s0, 0x0028($sp)           
/* 0445C 80863AAC 27BD0050 */  addiu   $sp, $sp, 0x0050           ## $sp = 00000000
/* 04460 80863AB0 03E00008 */  jr      $ra                        
/* 04464 80863AB4 00000000 */  nop


