glabel EnTp_Draw
/* 019CC 80B227AC 27BDFFA0 */  addiu   $sp, $sp, 0xFFA0           ## $sp = FFFFFFA0
/* 019D0 80B227B0 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 019D4 80B227B4 AFA40060 */  sw      $a0, 0x0060($sp)           
/* 019D8 80B227B8 AFA50064 */  sw      $a1, 0x0064($sp)           
/* 019DC 80B227BC 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 019E0 80B227C0 3C0680B2 */  lui     $a2, %hi(D_80B22B40)       ## $a2 = 80B20000
/* 019E4 80B227C4 24C62B40 */  addiu   $a2, $a2, %lo(D_80B22B40)  ## $a2 = 80B22B40
/* 019E8 80B227C8 27A40044 */  addiu   $a0, $sp, 0x0044           ## $a0 = FFFFFFE4
/* 019EC 80B227CC 240705AB */  addiu   $a3, $zero, 0x05AB         ## $a3 = 000005AB
/* 019F0 80B227D0 0C031AB1 */  jal     func_800C6AC4              
/* 019F4 80B227D4 AFA50054 */  sw      $a1, 0x0054($sp)           
/* 019F8 80B227D8 8FAF0060 */  lw      $t7, 0x0060($sp)           
/* 019FC 80B227DC 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 01A00 80B227E0 8FA80054 */  lw      $t0, 0x0054($sp)           
/* 01A04 80B227E4 8DF80150 */  lw      $t8, 0x0150($t7)           ## 00000150
/* 01A08 80B227E8 5301008A */  beql    $t8, $at, .L80B22A14       
/* 01A0C 80B227EC 8FAC0064 */  lw      $t4, 0x0064($sp)           
/* 01A10 80B227F0 85E2001C */  lh      $v0, 0x001C($t7)           ## 0000001C
/* 01A14 80B227F4 2401000C */  addiu   $at, $zero, 0x000C         ## $at = 0000000C
/* 01A18 80B227F8 8FB90064 */  lw      $t9, 0x0064($sp)           
/* 01A1C 80B227FC 04420004 */  bltzl   $v0, .L80B22810            
/* 01A20 80B22800 8F240000 */  lw      $a0, 0x0000($t9)           ## 00000000
/* 01A24 80B22804 14410025 */  bne     $v0, $at, .L80B2289C       
/* 01A28 80B22808 8FB80064 */  lw      $t8, 0x0064($sp)           
/* 01A2C 80B2280C 8F240000 */  lw      $a0, 0x0000($t9)           ## 00000000
.L80B22810:
/* 01A30 80B22810 0C024F46 */  jal     func_80093D18              
/* 01A34 80B22814 AFA80054 */  sw      $t0, 0x0054($sp)           
/* 01A38 80B22818 8FA80054 */  lw      $t0, 0x0054($sp)           
/* 01A3C 80B2281C 3C0ADA38 */  lui     $t2, 0xDA38                ## $t2 = DA380000
/* 01A40 80B22820 354A0003 */  ori     $t2, $t2, 0x0003           ## $t2 = DA380003
/* 01A44 80B22824 8D0202C0 */  lw      $v0, 0x02C0($t0)           ## 000002C0
/* 01A48 80B22828 3C0580B2 */  lui     $a1, %hi(D_80B22B50)       ## $a1 = 80B20000
/* 01A4C 80B2282C 24A52B50 */  addiu   $a1, $a1, %lo(D_80B22B50)  ## $a1 = 80B22B50
/* 01A50 80B22830 24490008 */  addiu   $t1, $v0, 0x0008           ## $t1 = 00000008
/* 01A54 80B22834 AD0902C0 */  sw      $t1, 0x02C0($t0)           ## 000002C0
/* 01A58 80B22838 AC4A0000 */  sw      $t2, 0x0000($v0)           ## 00000000
/* 01A5C 80B2283C 8FAB0064 */  lw      $t3, 0x0064($sp)           
/* 01A60 80B22840 240605B3 */  addiu   $a2, $zero, 0x05B3         ## $a2 = 000005B3
/* 01A64 80B22844 8D640000 */  lw      $a0, 0x0000($t3)           ## 00000000
/* 01A68 80B22848 AFA80054 */  sw      $t0, 0x0054($sp)           
/* 01A6C 80B2284C 0C0346A2 */  jal     Matrix_NewMtx              
/* 01A70 80B22850 AFA20040 */  sw      $v0, 0x0040($sp)           
/* 01A74 80B22854 8FA30040 */  lw      $v1, 0x0040($sp)           
/* 01A78 80B22858 8FA80054 */  lw      $t0, 0x0054($sp)           
/* 01A7C 80B2285C 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
/* 01A80 80B22860 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 01A84 80B22864 8D0202C0 */  lw      $v0, 0x02C0($t0)           ## 000002C0
/* 01A88 80B22868 3C0E0600 */  lui     $t6, 0x0600                ## $t6 = 06000000
/* 01A8C 80B2286C 25CE08D0 */  addiu   $t6, $t6, 0x08D0           ## $t6 = 060008D0
/* 01A90 80B22870 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 01A94 80B22874 3C0DDE00 */  lui     $t5, 0xDE00                ## $t5 = DE000000
/* 01A98 80B22878 AD0C02C0 */  sw      $t4, 0x02C0($t0)           ## 000002C0
/* 01A9C 80B2287C 3C064100 */  lui     $a2, 0x4100                ## $a2 = 41000000
/* 01AA0 80B22880 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 01AA4 80B22884 46006386 */  mov.s   $f14, $f12                 
/* 01AA8 80B22888 AC4E0004 */  sw      $t6, 0x0004($v0)           ## 00000004
/* 01AAC 80B2288C 0C034261 */  jal     Matrix_Translate              
/* 01AB0 80B22890 AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
/* 01AB4 80B22894 1000005F */  beq     $zero, $zero, .L80B22A14   
/* 01AB8 80B22898 8FAC0064 */  lw      $t4, 0x0064($sp)           
.L80B2289C:
/* 01ABC 80B2289C 8F040000 */  lw      $a0, 0x0000($t8)           ## 00000000
/* 01AC0 80B228A0 0C024F61 */  jal     func_80093D84              
/* 01AC4 80B228A4 AFA80054 */  sw      $t0, 0x0054($sp)           
/* 01AC8 80B228A8 8FA40064 */  lw      $a0, 0x0064($sp)           
/* 01ACC 80B228AC 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 01AD0 80B228B0 34211DA0 */  ori     $at, $at, 0x1DA0           ## $at = 00011DA0
/* 01AD4 80B228B4 0C0347F5 */  jal     func_800D1FD4              
/* 01AD8 80B228B8 00812021 */  addu    $a0, $a0, $at              
/* 01ADC 80B228BC 8FA80054 */  lw      $t0, 0x0054($sp)           
/* 01AE0 80B228C0 3C19FA00 */  lui     $t9, 0xFA00                ## $t9 = FA000000
/* 01AE4 80B228C4 3C05E700 */  lui     $a1, 0xE700                ## $a1 = E7000000
/* 01AE8 80B228C8 8D0202D0 */  lw      $v0, 0x02D0($t0)           ## 000002D0
/* 01AEC 80B228CC 3C040600 */  lui     $a0, 0x0600                ## $a0 = 06000000
/* 01AF0 80B228D0 24840C68 */  addiu   $a0, $a0, 0x0C68           ## $a0 = 06000C68
/* 01AF4 80B228D4 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 01AF8 80B228D8 AD0F02D0 */  sw      $t7, 0x02D0($t0)           ## 000002D0
/* 01AFC 80B228DC AC590000 */  sw      $t9, 0x0000($v0)           ## 00000000
/* 01B00 80B228E0 8FA90060 */  lw      $t1, 0x0060($sp)           
/* 01B04 80B228E4 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 01B08 80B228E8 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 01B0C 80B228EC 852B0160 */  lh      $t3, 0x0160($t1)           ## 00000160
/* 01B10 80B228F0 852E015E */  lh      $t6, 0x015E($t1)           ## 0000015E
/* 01B14 80B228F4 240605C8 */  addiu   $a2, $zero, 0x05C8         ## $a2 = 000005C8
/* 01B18 80B228F8 000B6600 */  sll     $t4, $t3, 24               
/* 01B1C 80B228FC 358DFF00 */  ori     $t5, $t4, 0xFF00           ## $t5 = 0000FF00
/* 01B20 80B22900 31D800FF */  andi    $t8, $t6, 0x00FF           ## $t8 = 00000000
/* 01B24 80B22904 01B87825 */  or      $t7, $t5, $t8              ## $t7 = 0000FF00
/* 01B28 80B22908 AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000004
/* 01B2C 80B2290C 8D0202D0 */  lw      $v0, 0x02D0($t0)           ## 000002D0
/* 01B30 80B22910 3C0C5566 */  lui     $t4, 0x5566                ## $t4 = 55660000
/* 01B34 80B22914 3C0BFC30 */  lui     $t3, 0xFC30                ## $t3 = FC300000
/* 01B38 80B22918 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 01B3C 80B2291C AD1902D0 */  sw      $t9, 0x02D0($t0)           ## 000002D0
/* 01B40 80B22920 AC450000 */  sw      $a1, 0x0000($v0)           ## 00000000
/* 01B44 80B22924 AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 01B48 80B22928 8D0202D0 */  lw      $v0, 0x02D0($t0)           ## 000002D0
/* 01B4C 80B2292C 356BB261 */  ori     $t3, $t3, 0xB261           ## $t3 = FC30B261
/* 01B50 80B22930 358CDB6D */  ori     $t4, $t4, 0xDB6D           ## $t4 = 5566DB6D
/* 01B54 80B22934 244A0008 */  addiu   $t2, $v0, 0x0008           ## $t2 = 00000008
/* 01B58 80B22938 AD0A02D0 */  sw      $t2, 0x02D0($t0)           ## 000002D0
/* 01B5C 80B2293C AC4B0000 */  sw      $t3, 0x0000($v0)           ## 00000000
/* 01B60 80B22940 AC4C0004 */  sw      $t4, 0x0004($v0)           ## 00000004
/* 01B64 80B22944 8D0202D0 */  lw      $v0, 0x02D0($t0)           ## 000002D0
/* 01B68 80B22948 0004C100 */  sll     $t8, $a0,  4               
/* 01B6C 80B2294C 00187F02 */  srl     $t7, $t8, 28               
/* 01B70 80B22950 24490008 */  addiu   $t1, $v0, 0x0008           ## $t1 = 00000008
/* 01B74 80B22954 AD0902D0 */  sw      $t1, 0x02D0($t0)           ## 000002D0
/* 01B78 80B22958 AC450000 */  sw      $a1, 0x0000($v0)           ## 00000000
/* 01B7C 80B2295C AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 01B80 80B22960 8D0202D0 */  lw      $v0, 0x02D0($t0)           ## 000002D0
/* 01B84 80B22964 3C0DDB06 */  lui     $t5, 0xDB06                ## $t5 = DB060000
/* 01B88 80B22968 35AD0020 */  ori     $t5, $t5, 0x0020           ## $t5 = DB060020
/* 01B8C 80B2296C 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 01B90 80B22970 000FC880 */  sll     $t9, $t7,  2               
/* 01B94 80B22974 3C0A8016 */  lui     $t2, 0x8016                ## $t2 = 80160000
/* 01B98 80B22978 AD0E02D0 */  sw      $t6, 0x02D0($t0)           ## 000002D0
/* 01B9C 80B2297C 01595021 */  addu    $t2, $t2, $t9              
/* 01BA0 80B22980 AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
/* 01BA4 80B22984 8D4A6FA8 */  lw      $t2, 0x6FA8($t2)           ## 80166FA8
/* 01BA8 80B22988 00815824 */  and     $t3, $a0, $at              
/* 01BAC 80B2298C 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 01BB0 80B22990 014B6021 */  addu    $t4, $t2, $t3              
/* 01BB4 80B22994 01814821 */  addu    $t1, $t4, $at              
/* 01BB8 80B22998 AC490004 */  sw      $t1, 0x0004($v0)           ## 00000004
/* 01BBC 80B2299C 8D0202D0 */  lw      $v0, 0x02D0($t0)           ## 000002D0
/* 01BC0 80B229A0 3C18DA38 */  lui     $t8, 0xDA38                ## $t8 = DA380000
/* 01BC4 80B229A4 37180003 */  ori     $t8, $t8, 0x0003           ## $t8 = DA380003
/* 01BC8 80B229A8 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 01BCC 80B229AC AD0E02D0 */  sw      $t6, 0x02D0($t0)           ## 000002D0
/* 01BD0 80B229B0 AC450000 */  sw      $a1, 0x0000($v0)           ## 00000000
/* 01BD4 80B229B4 AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 01BD8 80B229B8 8D0202D0 */  lw      $v0, 0x02D0($t0)           ## 000002D0
/* 01BDC 80B229BC 3C0580B2 */  lui     $a1, %hi(D_80B22B60)       ## $a1 = 80B20000
/* 01BE0 80B229C0 24A52B60 */  addiu   $a1, $a1, %lo(D_80B22B60)  ## $a1 = 80B22B60
/* 01BE4 80B229C4 244D0008 */  addiu   $t5, $v0, 0x0008           ## $t5 = 00000008
/* 01BE8 80B229C8 AD0D02D0 */  sw      $t5, 0x02D0($t0)           ## 000002D0
/* 01BEC 80B229CC AC580000 */  sw      $t8, 0x0000($v0)           ## 00000000
/* 01BF0 80B229D0 8FAF0064 */  lw      $t7, 0x0064($sp)           
/* 01BF4 80B229D4 8DE40000 */  lw      $a0, 0x0000($t7)           ## 0000FF00
/* 01BF8 80B229D8 AFA80054 */  sw      $t0, 0x0054($sp)           
/* 01BFC 80B229DC 0C0346A2 */  jal     Matrix_NewMtx              
/* 01C00 80B229E0 AFA20020 */  sw      $v0, 0x0020($sp)           
/* 01C04 80B229E4 8FA30020 */  lw      $v1, 0x0020($sp)           
/* 01C08 80B229E8 8FA80054 */  lw      $t0, 0x0054($sp)           
/* 01C0C 80B229EC 3C0B0600 */  lui     $t3, 0x0600                ## $t3 = 06000000
/* 01C10 80B229F0 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 01C14 80B229F4 8D0202D0 */  lw      $v0, 0x02D0($t0)           ## 000002D0
/* 01C18 80B229F8 256B0000 */  addiu   $t3, $t3, 0x0000           ## $t3 = 06000000
/* 01C1C 80B229FC 3C0ADE00 */  lui     $t2, 0xDE00                ## $t2 = DE000000
/* 01C20 80B22A00 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 01C24 80B22A04 AD1902D0 */  sw      $t9, 0x02D0($t0)           ## 000002D0
/* 01C28 80B22A08 AC4B0004 */  sw      $t3, 0x0004($v0)           ## 00000004
/* 01C2C 80B22A0C AC4A0000 */  sw      $t2, 0x0000($v0)           ## 00000000
/* 01C30 80B22A10 8FAC0064 */  lw      $t4, 0x0064($sp)           
.L80B22A14:
/* 01C34 80B22A14 3C0680B2 */  lui     $a2, %hi(D_80B22B70)       ## $a2 = 80B20000
/* 01C38 80B22A18 24C62B70 */  addiu   $a2, $a2, %lo(D_80B22B70)  ## $a2 = 80B22B70
/* 01C3C 80B22A1C 27A40044 */  addiu   $a0, $sp, 0x0044           ## $a0 = FFFFFFE4
/* 01C40 80B22A20 240705D7 */  addiu   $a3, $zero, 0x05D7         ## $a3 = 000005D7
/* 01C44 80B22A24 0C031AD5 */  jal     func_800C6B54              
/* 01C48 80B22A28 8D850000 */  lw      $a1, 0x0000($t4)           ## 00000000
/* 01C4C 80B22A2C 8FA30060 */  lw      $v1, 0x0060($sp)           
/* 01C50 80B22A30 2401000B */  addiu   $at, $zero, 0x000B         ## $at = 0000000B
/* 01C54 80B22A34 00002025 */  or      $a0, $zero, $zero          ## $a0 = 00000000
/* 01C58 80B22A38 8462001C */  lh      $v0, 0x001C($v1)           ## 0000001C
/* 01C5C 80B22A3C 18400003 */  blez    $v0, .L80B22A4C            
/* 01C60 80B22A40 00000000 */  nop
/* 01C64 80B22A44 54410004 */  bnel    $v0, $at, .L80B22A58       
/* 01C68 80B22A48 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80B22A4C:
/* 01C6C 80B22A4C 0C018A29 */  jal     func_800628A4              
/* 01C70 80B22A50 24650174 */  addiu   $a1, $v1, 0x0174           ## $a1 = 00000174
/* 01C74 80B22A54 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80B22A58:
/* 01C78 80B22A58 27BD0060 */  addiu   $sp, $sp, 0x0060           ## $sp = 00000000
/* 01C7C 80B22A5C 03E00008 */  jr      $ra                        
/* 01C80 80B22A60 00000000 */  nop
/* 01C84 80B22A64 00000000 */  nop
/* 01C88 80B22A68 00000000 */  nop
/* 01C8C 80B22A6C 00000000 */  nop

