glabel func_80997220
/* 00F80 80997220 27BDFFB0 */  addiu   $sp, $sp, 0xFFB0           ## $sp = FFFFFFB0
/* 00F84 80997224 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 00F88 80997228 AFB10018 */  sw      $s1, 0x0018($sp)           
/* 00F8C 8099722C AFB00014 */  sw      $s0, 0x0014($sp)           
/* 00F90 80997230 8CAE1C44 */  lw      $t6, 0x1C44($a1)           ## 00001C44
/* 00F94 80997234 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 00F98 80997238 00808825 */  or      $s1, $a0, $zero            ## $s1 = 00000000
/* 00F9C 8099723C AFAE004C */  sw      $t6, 0x004C($sp)           
/* 00FA0 80997240 80820003 */  lb      $v0, 0x0003($a0)           ## 00000003
/* 00FA4 80997244 27A5003C */  addiu   $a1, $sp, 0x003C           ## $a1 = FFFFFFEC
/* 00FA8 80997248 25C60024 */  addiu   $a2, $t6, 0x0024           ## $a2 = 00000024
/* 00FAC 8099724C 00021E00 */  sll     $v1, $v0, 24               
/* 00FB0 80997250 0440004C */  bltz    $v0, .L80997384            
/* 00FB4 80997254 00031E03 */  sra     $v1, $v1, 24               
/* 00FB8 80997258 0C00B6F4 */  jal     func_8002DBD0              
/* 00FBC 8099725C A3A3004B */  sb      $v1, 0x004B($sp)           
/* 00FC0 80997260 C7A40044 */  lwc1    $f4, 0x0044($sp)           
/* 00FC4 80997264 44803000 */  mtc1    $zero, $f6                 ## $f6 = 0.00
/* 00FC8 80997268 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00FCC 8099726C 02012021 */  addu    $a0, $s0, $at              
/* 00FD0 80997270 4606203C */  c.lt.s  $f4, $f6                   
/* 00FD4 80997274 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00FD8 80997278 83A3004B */  lb      $v1, 0x004B($sp)           
/* 00FDC 8099727C 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
/* 00FE0 80997280 45000003 */  bc1f    .L80997290                 
/* 00FE4 80997284 34211CBC */  ori     $at, $at, 0x1CBC           ## $at = 00011CBC
/* 00FE8 80997288 10000001 */  beq     $zero, $zero, .L80997290   
/* 00FEC 8099728C 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L80997290:
/* 00FF0 80997290 9638001C */  lhu     $t8, 0x001C($s1)           ## 0000001C
/* 00FF4 80997294 8C8F1D38 */  lw      $t7, 0x1D38($a0)           ## 00001D38
/* 00FF8 80997298 00025040 */  sll     $t2, $v0,  1               
/* 00FFC 8099729C 0018CA83 */  sra     $t9, $t8, 10               
/* 01000 809972A0 00194100 */  sll     $t0, $t9,  4               
/* 01004 809972A4 01E84821 */  addu    $t1, $t7, $t0              
/* 01008 809972A8 012A5821 */  addu    $t3, $t1, $t2              
/* 0100C 809972AC 816C0000 */  lb      $t4, 0x0000($t3)           ## 00000000
/* 01010 809972B0 02011021 */  addu    $v0, $s0, $at              
/* 01014 809972B4 A22C0003 */  sb      $t4, 0x0003($s1)           ## 00000003
/* 01018 809972B8 822D0003 */  lb      $t5, 0x0003($s1)           ## 00000003
/* 0101C 809972BC 506D0029 */  beql    $v1, $t5, .L80997364       
/* 01020 809972C0 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 01024 809972C4 8C580000 */  lw      $t8, 0x0000($v0)           ## 00000000
/* 01028 809972C8 27A50028 */  addiu   $a1, $sp, 0x0028           ## $a1 = FFFFFFD8
/* 0102C 809972CC 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 01030 809972D0 ACB80000 */  sw      $t8, 0x0000($a1)           ## FFFFFFD8
/* 01034 809972D4 8C4E0004 */  lw      $t6, 0x0004($v0)           ## 00000004
/* 01038 809972D8 34211CD0 */  ori     $at, $at, 0x1CD0           ## $at = 00011CD0
/* 0103C 809972DC 02011821 */  addu    $v1, $s0, $at              
/* 01040 809972E0 ACAE0004 */  sw      $t6, 0x0004($a1)           ## FFFFFFDC
/* 01044 809972E4 8C580008 */  lw      $t8, 0x0008($v0)           ## 00000008
/* 01048 809972E8 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 0104C 809972EC 00300821 */  addu    $at, $at, $s0              
/* 01050 809972F0 ACB80008 */  sw      $t8, 0x0008($a1)           ## FFFFFFE0
/* 01054 809972F4 8C4E000C */  lw      $t6, 0x000C($v0)           ## 0000000C
/* 01058 809972F8 ACAE000C */  sw      $t6, 0x000C($a1)           ## FFFFFFE4
/* 0105C 809972FC 8C580010 */  lw      $t8, 0x0010($v0)           ## 00000010
/* 01060 80997300 ACB80010 */  sw      $t8, 0x0010($a1)           ## FFFFFFE8
/* 01064 80997304 8C6F0000 */  lw      $t7, 0x0000($v1)           ## 00000000
/* 01068 80997308 AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
/* 0106C 8099730C 8C790004 */  lw      $t9, 0x0004($v1)           ## 00000004
/* 01070 80997310 AC590004 */  sw      $t9, 0x0004($v0)           ## 00000004
/* 01074 80997314 8C6F0008 */  lw      $t7, 0x0008($v1)           ## 00000008
/* 01078 80997318 AC4F0008 */  sw      $t7, 0x0008($v0)           ## 00000008
/* 0107C 8099731C 8C79000C */  lw      $t9, 0x000C($v1)           ## 0000000C
/* 01080 80997320 AC59000C */  sw      $t9, 0x000C($v0)           ## 0000000C
/* 01084 80997324 8C6F0010 */  lw      $t7, 0x0010($v1)           ## 00000010
/* 01088 80997328 AC4F0010 */  sw      $t7, 0x0010($v0)           ## 00000010
/* 0108C 8099732C 8CA90000 */  lw      $t1, 0x0000($a1)           ## FFFFFFD8
/* 01090 80997330 AC690000 */  sw      $t1, 0x0000($v1)           ## 00000000
/* 01094 80997334 8CA80004 */  lw      $t0, 0x0004($a1)           ## FFFFFFDC
/* 01098 80997338 AC680004 */  sw      $t0, 0x0004($v1)           ## 00000004
/* 0109C 8099733C 8CA90008 */  lw      $t1, 0x0008($a1)           ## FFFFFFE0
/* 010A0 80997340 AC690008 */  sw      $t1, 0x0008($v1)           ## 00000008
/* 010A4 80997344 8CA8000C */  lw      $t0, 0x000C($a1)           ## FFFFFFE4
/* 010A8 80997348 AC68000C */  sw      $t0, 0x000C($v1)           ## 0000000C
/* 010AC 8099734C 8CA90010 */  lw      $t1, 0x0010($a1)           ## FFFFFFE8
/* 010B0 80997350 AC690010 */  sw      $t1, 0x0010($v1)           ## 00000010
/* 010B4 80997354 908A1CEC */  lbu     $t2, 0x1CEC($a0)           ## 00001CEC
/* 010B8 80997358 394B0001 */  xori    $t3, $t2, 0x0001           ## $t3 = 00000001
/* 010BC 8099735C A02B1CEC */  sb      $t3, 0x1CEC($at)           ## 00011CEC
/* 010C0 80997360 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
.L80997364:
/* 010C4 80997364 34211CBC */  ori     $at, $at, 0x1CBC           ## $at = 00011CBC
/* 010C8 80997368 02012821 */  addu    $a1, $s0, $at              
/* 010CC 8099736C 0C025D4D */  jal     func_80097534              
/* 010D0 80997370 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 010D4 80997374 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 010D8 80997378 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 010DC 8099737C 0C0302BD */  jal     Gameplay_SetupRespawnPoint              
/* 010E0 80997380 24060EFF */  addiu   $a2, $zero, 0x0EFF         ## $a2 = 00000EFF
.L80997384:
/* 010E4 80997384 44804000 */  mtc1    $zero, $f8                 ## $f8 = 0.00
/* 010E8 80997388 A6200164 */  sh      $zero, 0x0164($s1)         ## 00000164
/* 010EC 8099738C 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 010F0 80997390 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 010F4 80997394 0C2658AB */  jal     func_809962AC              
/* 010F8 80997398 E6280060 */  swc1    $f8, 0x0060($s1)           ## 00000060
/* 010FC 8099739C 1040000D */  beq     $v0, $zero, .L809973D4     
/* 01100 809973A0 8FAC004C */  lw      $t4, 0x004C($sp)           
/* 01104 809973A4 8D8D067C */  lw      $t5, 0x067C($t4)           ## 0000067C
/* 01108 809973A8 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 0110C 809973AC 3C058099 */  lui     $a1, %hi(func_80997568)    ## $a1 = 80990000
/* 01110 809973B0 31AE0800 */  andi    $t6, $t5, 0x0800           ## $t6 = 00000000
/* 01114 809973B4 55C00008 */  bnel    $t6, $zero, .L809973D8     
/* 01118 809973B8 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 0111C 809973BC 0C2658A8 */  jal     func_809962A0              
/* 01120 809973C0 24A57568 */  addiu   $a1, $a1, %lo(func_80997568) ## $a1 = 80997568
/* 01124 809973C4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01128 809973C8 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 0112C 809973CC 0C00B7D5 */  jal     func_8002DF54              
/* 01130 809973D0 24060002 */  addiu   $a2, $zero, 0x0002         ## $a2 = 00000002
.L809973D4:
/* 01134 809973D4 8FBF001C */  lw      $ra, 0x001C($sp)           
.L809973D8:
/* 01138 809973D8 8FB00014 */  lw      $s0, 0x0014($sp)           
/* 0113C 809973DC 8FB10018 */  lw      $s1, 0x0018($sp)           
/* 01140 809973E0 03E00008 */  jr      $ra                        
/* 01144 809973E4 27BD0050 */  addiu   $sp, $sp, 0x0050           ## $sp = 00000000
