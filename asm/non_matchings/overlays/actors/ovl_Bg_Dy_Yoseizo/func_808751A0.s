glabel func_808751A0
/* 02970 808751A0 27BDFF60 */  addiu   $sp, $sp, 0xFF60           ## $sp = FFFFFF60
/* 02974 808751A4 AFBF0044 */  sw      $ra, 0x0044($sp)           
/* 02978 808751A8 AFBE0040 */  sw      $s8, 0x0040($sp)           
/* 0297C 808751AC AFB7003C */  sw      $s7, 0x003C($sp)           
/* 02980 808751B0 AFB60038 */  sw      $s6, 0x0038($sp)           
/* 02984 808751B4 AFB50034 */  sw      $s5, 0x0034($sp)           
/* 02988 808751B8 AFB40030 */  sw      $s4, 0x0030($sp)           
/* 0298C 808751BC AFB3002C */  sw      $s3, 0x002C($sp)           
/* 02990 808751C0 AFB20028 */  sw      $s2, 0x0028($sp)           
/* 02994 808751C4 AFB10024 */  sw      $s1, 0x0024($sp)           
/* 02998 808751C8 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 0299C 808751CC F7B40018 */  sdc1    $f20, 0x0018($sp)          
/* 029A0 808751D0 AFA500A4 */  sw      $a1, 0x00A4($sp)           
/* 029A4 808751D4 8CB10000 */  lw      $s1, 0x0000($a1)           ## 00000000
/* 029A8 808751D8 24900394 */  addiu   $s0, $a0, 0x0394           ## $s0 = 00000394
/* 029AC 808751DC 3C068087 */  lui     $a2, %hi(D_80875638)       ## $a2 = 80870000
/* 029B0 808751E0 00009825 */  or      $s3, $zero, $zero          ## $s3 = 00000000
/* 029B4 808751E4 24C65638 */  addiu   $a2, $a2, %lo(D_80875638)  ## $a2 = 80875638
/* 029B8 808751E8 27A4007C */  addiu   $a0, $sp, 0x007C           ## $a0 = FFFFFFDC
/* 029BC 808751EC 240706E7 */  addiu   $a3, $zero, 0x06E7         ## $a3 = 000006E7
/* 029C0 808751F0 0C031AB1 */  jal     func_800C6AC4              
/* 029C4 808751F4 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 029C8 808751F8 8FAF00A4 */  lw      $t7, 0x00A4($sp)           
/* 029CC 808751FC 0C024F61 */  jal     func_80093D84              
/* 029D0 80875200 8DE40000 */  lw      $a0, 0x0000($t7)           ## 00000000
/* 029D4 80875204 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 029D8 80875208 4481A000 */  mtc1    $at, $f20                  ## $f20 = 1.00
/* 029DC 8087520C 0000A025 */  or      $s4, $zero, $zero          ## $s4 = 00000000
/* 029E0 80875210 3C1E8000 */  lui     $s8, 0x8000                ## $s8 = 80000000
/* 029E4 80875214 3C17DE00 */  lui     $s7, 0xDE00                ## $s7 = DE000000
.L80875218:
/* 029E8 80875218 92180000 */  lbu     $t8, 0x0000($s0)           ## 00000394
/* 029EC 8087521C 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 029F0 80875220 8FB200A4 */  lw      $s2, 0x00A4($sp)           
/* 029F4 80875224 1701006B */  bne     $t8, $at, .L808753D4       
/* 029F8 80875228 3C020600 */  lui     $v0, 0x0600                ## $v0 = 06000000
/* 029FC 8087522C 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 02A00 80875230 34211DA0 */  ori     $at, $at, 0x1DA0           ## $at = 00011DA0
/* 02A04 80875234 244258D8 */  addiu   $v0, $v0, 0x58D8           ## $v0 = 060058D8
/* 02A08 80875238 0002C900 */  sll     $t9, $v0,  4               
/* 02A0C 8087523C 02419021 */  addu    $s2, $s2, $at              
/* 02A10 80875240 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 02A14 80875244 00194702 */  srl     $t0, $t9, 28               
/* 02A18 80875248 3C0A8016 */  lui     $t2, 0x8016                ## $t2 = 80160000
/* 02A1C 8087524C 254A6FA8 */  addiu   $t2, $t2, 0x6FA8           ## $t2 = 80166FA8
/* 02A20 80875250 00084880 */  sll     $t1, $t0,  2               
/* 02A24 80875254 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 02A28 80875258 0041B024 */  and     $s6, $v0, $at              
/* 02A2C 8087525C 1660001A */  bne     $s3, $zero, .L808752C8     
/* 02A30 80875260 012AA821 */  addu    $s5, $t1, $t2              
/* 02A34 80875264 3C020600 */  lui     $v0, 0x0600                ## $v0 = 06000000
/* 02A38 80875268 24425860 */  addiu   $v0, $v0, 0x5860           ## $v0 = 06005860
/* 02A3C 8087526C 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 02A40 80875270 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 02A44 80875274 00412824 */  and     $a1, $v0, $at              
/* 02A48 80875278 00025900 */  sll     $t3, $v0,  4               
/* 02A4C 8087527C 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 02A50 80875280 000B6702 */  srl     $t4, $t3, 28               
/* 02A54 80875284 000C6880 */  sll     $t5, $t4,  2               
/* 02A58 80875288 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 06005868
/* 02A5C 8087528C AE2E02D0 */  sw      $t6, 0x02D0($s1)           ## 000002D0
/* 02A60 80875290 01AA2021 */  addu    $a0, $t5, $t2              
/* 02A64 80875294 AC570000 */  sw      $s7, 0x0000($v0)           ## 06005860
/* 02A68 80875298 8C8F0000 */  lw      $t7, 0x0000($a0)           ## 00000000
/* 02A6C 8087529C 3C09E700 */  lui     $t1, 0xE700                ## $t1 = E7000000
/* 02A70 808752A0 26730001 */  addiu   $s3, $s3, 0x0001           ## $s3 = 00000001
/* 02A74 808752A4 01E5C021 */  addu    $t8, $t7, $a1              
/* 02A78 808752A8 031EC821 */  addu    $t9, $t8, $s8              
/* 02A7C 808752AC AC590004 */  sw      $t9, 0x0004($v0)           ## 06005864
/* 02A80 808752B0 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 02A84 808752B4 327300FF */  andi    $s3, $s3, 0x00FF           ## $s3 = 00000001
/* 02A88 808752B8 24480008 */  addiu   $t0, $v0, 0x0008           ## $t0 = 06005868
/* 02A8C 808752BC AE2802D0 */  sw      $t0, 0x02D0($s1)           ## 000002D0
/* 02A90 808752C0 AC400004 */  sw      $zero, 0x0004($v0)         ## 06005864
/* 02A94 808752C4 AC490000 */  sw      $t1, 0x0000($v0)           ## 06005860
.L808752C8:
/* 02A98 808752C8 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 02A9C 808752CC 3C0CFA00 */  lui     $t4, 0xFA00                ## $t4 = FA000000
/* 02AA0 808752D0 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 02AA4 808752D4 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 06005868
/* 02AA8 808752D8 AE2B02D0 */  sw      $t3, 0x02D0($s1)           ## 000002D0
/* 02AAC 808752DC AC4C0000 */  sw      $t4, 0x0000($v0)           ## 06005860
/* 02AB0 808752E0 920A0028 */  lbu     $t2, 0x0028($s0)           ## 000003BC
/* 02AB4 808752E4 92180029 */  lbu     $t8, 0x0029($s0)           ## 000003BD
/* 02AB8 808752E8 920B002A */  lbu     $t3, 0x002A($s0)           ## 000003BE
/* 02ABC 808752EC 000A7600 */  sll     $t6, $t2, 24               
/* 02AC0 808752F0 860A002E */  lh      $t2, 0x002E($s0)           ## 000003C2
/* 02AC4 808752F4 0018CC00 */  sll     $t9, $t8, 16               
/* 02AC8 808752F8 01D94025 */  or      $t0, $t6, $t9              ## $t0 = 06005868
/* 02ACC 808752FC 000B6200 */  sll     $t4, $t3,  8               
/* 02AD0 80875300 010C6825 */  or      $t5, $t0, $t4              ## $t5 = FE005868
/* 02AD4 80875304 314F00FF */  andi    $t7, $t2, 0x00FF           ## $t7 = 000000A8
/* 02AD8 80875308 01AFC025 */  or      $t8, $t5, $t7              ## $t8 = FE0058E8
/* 02ADC 8087530C AC580004 */  sw      $t8, 0x0004($v0)           ## 06005864
/* 02AE0 80875310 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 02AE4 80875314 3C19FB00 */  lui     $t9, 0xFB00                ## $t9 = FB000000
/* 02AE8 80875318 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 06005868
/* 02AEC 8087531C AE2E02D0 */  sw      $t6, 0x02D0($s1)           ## 000002D0
/* 02AF0 80875320 AC590000 */  sw      $t9, 0x0000($v0)           ## 06005860
/* 02AF4 80875324 920A002C */  lbu     $t2, 0x002C($s0)           ## 000003C0
/* 02AF8 80875328 920B002B */  lbu     $t3, 0x002B($s0)           ## 000003BF
/* 02AFC 8087532C 920E002D */  lbu     $t6, 0x002D($s0)           ## 000003C1
/* 02B00 80875330 000A6C00 */  sll     $t5, $t2, 16               
/* 02B04 80875334 000B4600 */  sll     $t0, $t3, 24               
/* 02B08 80875338 010D7825 */  or      $t7, $t0, $t5              ## $t7 = FE005868
/* 02B0C 8087533C 000ECA00 */  sll     $t9, $t6,  8               
/* 02B10 80875340 01F94825 */  or      $t1, $t7, $t9              ## $t1 = FF005868
/* 02B14 80875344 AC490004 */  sw      $t1, 0x0004($v0)           ## 06005864
/* 02B18 80875348 8E06000C */  lw      $a2, 0x000C($s0)           ## 000003A0
/* 02B1C 8087534C C60E0008 */  lwc1    $f14, 0x0008($s0)          ## 0000039C
/* 02B20 80875350 0C034261 */  jal     Matrix_Translate              
/* 02B24 80875354 C60C0004 */  lwc1    $f12, 0x0004($s0)          ## 00000398
/* 02B28 80875358 0C0347F5 */  jal     func_800D1FD4              
/* 02B2C 8087535C 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 02B30 80875360 C60C0030 */  lwc1    $f12, 0x0030($s0)          ## 000003C4
/* 02B34 80875364 4406A000 */  mfc1    $a2, $f20                  
/* 02B38 80875368 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 02B3C 8087536C 0C0342A3 */  jal     Matrix_Scale              
/* 02B40 80875370 46006386 */  mov.s   $f14, $f12                 
/* 02B44 80875374 C60C0040 */  lwc1    $f12, 0x0040($s0)          ## 000003D4
/* 02B48 80875378 0C0343B5 */  jal     Matrix_RotateZ              
/* 02B4C 8087537C 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 02B50 80875380 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 02B54 80875384 3C0CDA38 */  lui     $t4, 0xDA38                ## $t4 = DA380000
/* 02B58 80875388 358C0003 */  ori     $t4, $t4, 0x0003           ## $t4 = DA380003
/* 02B5C 8087538C 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 02B60 80875390 AE2B02D0 */  sw      $t3, 0x02D0($s1)           ## 000002D0
/* 02B64 80875394 3C058087 */  lui     $a1, %hi(D_80875650)       ## $a1 = 80870000
/* 02B68 80875398 24A55650 */  addiu   $a1, $a1, %lo(D_80875650)  ## $a1 = 80875650
/* 02B6C 8087539C 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 02B70 808753A0 24060712 */  addiu   $a2, $zero, 0x0712         ## $a2 = 00000712
/* 02B74 808753A4 AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 02B78 808753A8 0C0346A2 */  jal     Matrix_NewMtx              
/* 02B7C 808753AC 00409025 */  or      $s2, $v0, $zero            ## $s2 = 00000000
/* 02B80 808753B0 AE420004 */  sw      $v0, 0x0004($s2)           ## 00000004
/* 02B84 808753B4 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 02B88 808753B8 244A0008 */  addiu   $t2, $v0, 0x0008           ## $t2 = 00000008
/* 02B8C 808753BC AE2A02D0 */  sw      $t2, 0x02D0($s1)           ## 000002D0
/* 02B90 808753C0 AC570000 */  sw      $s7, 0x0000($v0)           ## 00000000
/* 02B94 808753C4 8EA80000 */  lw      $t0, 0x0000($s5)           ## 00000000
/* 02B98 808753C8 01166821 */  addu    $t5, $t0, $s6              
/* 02B9C 808753CC 01BEC021 */  addu    $t8, $t5, $s8              
/* 02BA0 808753D0 AC580004 */  sw      $t8, 0x0004($v0)           ## 00000004
.L808753D4:
/* 02BA4 808753D4 26940001 */  addiu   $s4, $s4, 0x0001           ## $s4 = 00000001
/* 02BA8 808753D8 0014A400 */  sll     $s4, $s4, 16               
/* 02BAC 808753DC 0014A403 */  sra     $s4, $s4, 16               
/* 02BB0 808753E0 2A8100C8 */  slti    $at, $s4, 0x00C8           
/* 02BB4 808753E4 1420FF8C */  bne     $at, $zero, .L80875218     
/* 02BB8 808753E8 26100044 */  addiu   $s0, $s0, 0x0044           ## $s0 = 000003D8
/* 02BBC 808753EC 3C068087 */  lui     $a2, %hi(D_80875668)       ## $a2 = 80870000
/* 02BC0 808753F0 24C65668 */  addiu   $a2, $a2, %lo(D_80875668)  ## $a2 = 80875668
/* 02BC4 808753F4 27A4007C */  addiu   $a0, $sp, 0x007C           ## $a0 = FFFFFFDC
/* 02BC8 808753F8 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 02BCC 808753FC 0C031AD5 */  jal     func_800C6B54              
/* 02BD0 80875400 2407071B */  addiu   $a3, $zero, 0x071B         ## $a3 = 0000071B
/* 02BD4 80875404 8FBF0044 */  lw      $ra, 0x0044($sp)           
/* 02BD8 80875408 D7B40018 */  ldc1    $f20, 0x0018($sp)          
/* 02BDC 8087540C 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 02BE0 80875410 8FB10024 */  lw      $s1, 0x0024($sp)           
/* 02BE4 80875414 8FB20028 */  lw      $s2, 0x0028($sp)           
/* 02BE8 80875418 8FB3002C */  lw      $s3, 0x002C($sp)           
/* 02BEC 8087541C 8FB40030 */  lw      $s4, 0x0030($sp)           
/* 02BF0 80875420 8FB50034 */  lw      $s5, 0x0034($sp)           
/* 02BF4 80875424 8FB60038 */  lw      $s6, 0x0038($sp)           
/* 02BF8 80875428 8FB7003C */  lw      $s7, 0x003C($sp)           
/* 02BFC 8087542C 8FBE0040 */  lw      $s8, 0x0040($sp)           
/* 02C00 80875430 03E00008 */  jr      $ra                        
/* 02C04 80875434 27BD00A0 */  addiu   $sp, $sp, 0x00A0           ## $sp = 00000000
/* 02C08 80875438 00000000 */  nop
/* 02C0C 8087543C 00000000 */  nop

