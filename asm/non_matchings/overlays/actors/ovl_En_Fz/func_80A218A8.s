.rdata
glabel D_80A21C90
    .asciz "../z_en_fz.c"
    .balign 4

glabel D_80A21CA0
    .asciz "../z_en_fz.c"
    .balign 4

glabel D_80A21CB0
    .asciz "../z_en_fz.c"
    .balign 4

.text
glabel func_80A218A8
/* 01AE8 80A218A8 27BDFF40 */  addiu   $sp, $sp, 0xFF40           ## $sp = FFFFFF40
/* 01AEC 80A218AC AFBF0064 */  sw      $ra, 0x0064($sp)           
/* 01AF0 80A218B0 AFBE0060 */  sw      $s8, 0x0060($sp)           
/* 01AF4 80A218B4 AFB7005C */  sw      $s7, 0x005C($sp)           
/* 01AF8 80A218B8 AFB60058 */  sw      $s6, 0x0058($sp)           
/* 01AFC 80A218BC AFB50054 */  sw      $s5, 0x0054($sp)           
/* 01B00 80A218C0 AFB40050 */  sw      $s4, 0x0050($sp)           
/* 01B04 80A218C4 AFB3004C */  sw      $s3, 0x004C($sp)           
/* 01B08 80A218C8 AFB20048 */  sw      $s2, 0x0048($sp)           
/* 01B0C 80A218CC AFB10044 */  sw      $s1, 0x0044($sp)           
/* 01B10 80A218D0 AFB00040 */  sw      $s0, 0x0040($sp)           
/* 01B14 80A218D4 F7B40038 */  sdc1    $f20, 0x0038($sp)          
/* 01B18 80A218D8 AFA500C4 */  sw      $a1, 0x00C4($sp)           
/* 01B1C 80A218DC 8CB00000 */  lw      $s0, 0x0000($a1)           ## 00000000
/* 01B20 80A218E0 24920274 */  addiu   $s2, $a0, 0x0274           ## $s2 = 00000274
/* 01B24 80A218E4 3C0680A2 */  lui     $a2, %hi(D_80A21C90)       ## $a2 = 80A20000
/* 01B28 80A218E8 0000A825 */  or      $s5, $zero, $zero          ## $s5 = 00000000
/* 01B2C 80A218EC 24C61C90 */  addiu   $a2, $a2, %lo(D_80A21C90)  ## $a2 = 80A21C90
/* 01B30 80A218F0 27A4009C */  addiu   $a0, $sp, 0x009C           ## $a0 = FFFFFFDC
/* 01B34 80A218F4 24070568 */  addiu   $a3, $zero, 0x0568         ## $a3 = 00000568
/* 01B38 80A218F8 0C031AB1 */  jal     Graph_OpenDisps              
/* 01B3C 80A218FC 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 01B40 80A21900 8FAF00C4 */  lw      $t7, 0x00C4($sp)           
/* 01B44 80A21904 0C024F61 */  jal     func_80093D84              
/* 01B48 80A21908 8DE40000 */  lw      $a0, 0x0000($t7)           ## 00000000
/* 01B4C 80A2190C 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 01B50 80A21910 4481A000 */  mtc1    $at, $f20                  ## $f20 = 1.00
/* 01B54 80A21914 00009825 */  or      $s3, $zero, $zero          ## $s3 = 00000000
/* 01B58 80A21918 3C1EDE00 */  lui     $s8, 0xDE00                ## $s8 = DE000000
.L80A2191C:
/* 01B5C 80A2191C 92580000 */  lbu     $t8, 0x0000($s2)           ## 00000274
/* 01B60 80A21920 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 01B64 80A21924 34211DA0 */  ori     $at, $at, 0x1DA0           ## $at = 00011DA0
/* 01B68 80A21928 1B000078 */  blez    $t8, .L80A21B0C            
/* 01B6C 80A2192C 8FB400C4 */  lw      $s4, 0x00C4($sp)           
/* 01B70 80A21930 0281A021 */  addu    $s4, $s4, $at              
/* 01B74 80A21934 3C020600 */  lui     $v0, 0x0600                ## $v0 = 06000000
/* 01B78 80A21938 24423158 */  addiu   $v0, $v0, 0x3158           ## $v0 = 06003158
/* 01B7C 80A2193C 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 01B80 80A21940 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 01B84 80A21944 0041B824 */  and     $s7, $v0, $at              
/* 01B88 80A21948 0002C900 */  sll     $t9, $v0,  4               
/* 01B8C 80A2194C 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 01B90 80A21950 00194F02 */  srl     $t1, $t9, 28               
/* 01B94 80A21954 3C0B8016 */  lui     $t3, %hi(gSegments)
/* 01B98 80A21958 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 06003160
/* 01B9C 80A2195C AE0C02D0 */  sw      $t4, 0x02D0($s0)           ## 000002D0
/* 01BA0 80A21960 256B6FA8 */  addiu   $t3, %lo(gSegments)
/* 01BA4 80A21964 00095080 */  sll     $t2, $t1,  2               
/* 01BA8 80A21968 3C0DE700 */  lui     $t5, 0xE700                ## $t5 = E7000000
/* 01BAC 80A2196C 014BB021 */  addu    $s6, $t2, $t3              
/* 01BB0 80A21970 AC4D0000 */  sw      $t5, 0x0000($v0)           ## 06003158
/* 01BB4 80A21974 16A00017 */  bne     $s5, $zero, .L80A219D4     
/* 01BB8 80A21978 AC400004 */  sw      $zero, 0x0004($v0)         ## 0600315C
/* 01BBC 80A2197C 3C020600 */  lui     $v0, 0x0600                ## $v0 = 06000000
/* 01BC0 80A21980 244230A0 */  addiu   $v0, $v0, 0x30A0           ## $v0 = 060030A0
/* 01BC4 80A21984 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 01BC8 80A21988 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 01BCC 80A2198C 00412824 */  and     $a1, $v0, $at              
/* 01BD0 80A21990 00027100 */  sll     $t6, $v0,  4               
/* 01BD4 80A21994 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 01BD8 80A21998 000E7F02 */  srl     $t7, $t6, 28               
/* 01BDC 80A2199C 3C198016 */  lui     $t9, %hi(gSegments)
/* 01BE0 80A219A0 24490008 */  addiu   $t1, $v0, 0x0008           ## $t1 = 060030A8
/* 01BE4 80A219A4 AE0902D0 */  sw      $t1, 0x02D0($s0)           ## 000002D0
/* 01BE8 80A219A8 27396FA8 */  addiu   $t9, %lo(gSegments)
/* 01BEC 80A219AC 000FC080 */  sll     $t8, $t7,  2               
/* 01BF0 80A219B0 03192021 */  addu    $a0, $t8, $t9              
/* 01BF4 80A219B4 AC5E0000 */  sw      $s8, 0x0000($v0)           ## 060030A0
/* 01BF8 80A219B8 8C8A0000 */  lw      $t2, 0x0000($a0)           ## 00000000
/* 01BFC 80A219BC 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 01C00 80A219C0 26B50001 */  addiu   $s5, $s5, 0x0001           ## $s5 = 00000001
/* 01C04 80A219C4 01455821 */  addu    $t3, $t2, $a1              
/* 01C08 80A219C8 01616021 */  addu    $t4, $t3, $at              
/* 01C0C 80A219CC 32B500FF */  andi    $s5, $s5, 0x00FF           ## $s5 = 00000001
/* 01C10 80A219D0 AC4C0004 */  sw      $t4, 0x0004($v0)           ## 060030A4
.L80A219D4:
/* 01C14 80A219D4 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 01C18 80A219D8 3C0EFA00 */  lui     $t6, 0xFA00                ## $t6 = FA000000
/* 01C1C 80A219DC 3C01C3E1 */  lui     $at, 0xC3E1                ## $at = C3E10000
/* 01C20 80A219E0 244D0008 */  addiu   $t5, $v0, 0x0008           ## $t5 = 060030A8
/* 01C24 80A219E4 AE0D02D0 */  sw      $t5, 0x02D0($s0)           ## 000002D0
/* 01C28 80A219E8 AC4E0000 */  sw      $t6, 0x0000($v0)           ## 060030A0
/* 01C2C 80A219EC 864F002C */  lh      $t7, 0x002C($s2)           ## 000002A0
/* 01C30 80A219F0 24080003 */  addiu   $t0, $zero, 0x0003         ## $t0 = 00000003
/* 01C34 80A219F4 02680019 */  multu   $s3, $t0                   
/* 01C38 80A219F8 3421EB00 */  ori     $at, $at, 0xEB00           ## $at = C3E1EB00
/* 01C3C 80A219FC 31F800FF */  andi    $t8, $t7, 0x00FF           ## $t8 = 00000000
/* 01C40 80A21A00 0301C825 */  or      $t9, $t8, $at              ## $t9 = C3E1EB00
/* 01C44 80A21A04 AC590004 */  sw      $t9, 0x0004($v0)           ## 060030A4
/* 01C48 80A21A08 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 01C4C 80A21A0C 3C0ADB06 */  lui     $t2, 0xDB06                ## $t2 = DB060000
/* 01C50 80A21A10 354A0020 */  ori     $t2, $t2, 0x0020           ## $t2 = DB060020
/* 01C54 80A21A14 24490008 */  addiu   $t1, $v0, 0x0008           ## $t1 = 060030A8
/* 01C58 80A21A18 AE0902D0 */  sw      $t1, 0x02D0($s0)           ## 000002D0
/* 01C5C 80A21A1C AC4A0000 */  sw      $t2, 0x0000($v0)           ## 060030A0
/* 01C60 80A21A20 924C0001 */  lbu     $t4, 0x0001($s2)           ## 00000275
/* 01C64 80A21A24 00006812 */  mflo    $t5                        
/* 01C68 80A21A28 8FAB00C4 */  lw      $t3, 0x00C4($sp)           
/* 01C6C 80A21A2C 018D1821 */  addu    $v1, $t4, $t5              
/* 01C70 80A21A30 00680019 */  multu   $v1, $t0                   
/* 01C74 80A21A34 8D640000 */  lw      $a0, 0x0000($t3)           ## 80166FA8
/* 01C78 80A21A38 00033900 */  sll     $a3, $v1,  4               
/* 01C7C 80A21A3C 24090020 */  addiu   $t1, $zero, 0x0020         ## $t1 = 00000020
/* 01C80 80A21A40 24190020 */  addiu   $t9, $zero, 0x0020         ## $t9 = 00000020
/* 01C84 80A21A44 24180001 */  addiu   $t8, $zero, 0x0001         ## $t8 = 00000001
/* 01C88 80A21A48 240F0040 */  addiu   $t7, $zero, 0x0040         ## $t7 = 00000040
/* 01C8C 80A21A4C 240E0020 */  addiu   $t6, $zero, 0x0020         ## $t6 = 00000020
/* 01C90 80A21A50 AFAE0010 */  sw      $t6, 0x0010($sp)           
/* 01C94 80A21A54 AFAF0014 */  sw      $t7, 0x0014($sp)           
/* 01C98 80A21A58 00003012 */  mflo    $a2                        
/* 01C9C 80A21A5C AFB80018 */  sw      $t8, 0x0018($sp)           
/* 01CA0 80A21A60 AFB90024 */  sw      $t9, 0x0024($sp)           
/* 01CA4 80A21A64 AFA90028 */  sw      $t1, 0x0028($sp)           
/* 01CA8 80A21A68 00E33823 */  subu    $a3, $a3, $v1              
/* 01CAC 80A21A6C AFA00020 */  sw      $zero, 0x0020($sp)         
/* 01CB0 80A21A70 AFA0001C */  sw      $zero, 0x001C($sp)         
/* 01CB4 80A21A74 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 01CB8 80A21A78 0C0253D0 */  jal     Gfx_TwoTexScroll              
/* 01CBC 80A21A7C 00408825 */  or      $s1, $v0, $zero            ## $s1 = 060030A0
/* 01CC0 80A21A80 AE220004 */  sw      $v0, 0x0004($s1)           ## 060030A4
/* 01CC4 80A21A84 8E46000C */  lw      $a2, 0x000C($s2)           ## 00000280
/* 01CC8 80A21A88 C64E0008 */  lwc1    $f14, 0x0008($s2)          ## 0000027C
/* 01CCC 80A21A8C C64C0004 */  lwc1    $f12, 0x0004($s2)          ## 00000278
/* 01CD0 80A21A90 0C034261 */  jal     Matrix_Translate              
/* 01CD4 80A21A94 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 01CD8 80A21A98 0C0347F5 */  jal     func_800D1FD4              
/* 01CDC 80A21A9C 02802025 */  or      $a0, $s4, $zero            ## $a0 = 00000000
/* 01CE0 80A21AA0 C64C0030 */  lwc1    $f12, 0x0030($s2)          ## 000002A4
/* 01CE4 80A21AA4 4406A000 */  mfc1    $a2, $f20                  
/* 01CE8 80A21AA8 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 01CEC 80A21AAC 0C0342A3 */  jal     Matrix_Scale              
/* 01CF0 80A21AB0 46006386 */  mov.s   $f14, $f12                 
/* 01CF4 80A21AB4 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 01CF8 80A21AB8 3C0BDA38 */  lui     $t3, 0xDA38                ## $t3 = DA380000
/* 01CFC 80A21ABC 356B0003 */  ori     $t3, $t3, 0x0003           ## $t3 = DA380003
/* 01D00 80A21AC0 244A0008 */  addiu   $t2, $v0, 0x0008           ## $t2 = 00000008
/* 01D04 80A21AC4 AE0A02D0 */  sw      $t2, 0x02D0($s0)           ## 000002D0
/* 01D08 80A21AC8 3C0580A2 */  lui     $a1, %hi(D_80A21CA0)       ## $a1 = 80A20000
/* 01D0C 80A21ACC 24A51CA0 */  addiu   $a1, $a1, %lo(D_80A21CA0)  ## $a1 = 80A21CA0
/* 01D10 80A21AD0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01D14 80A21AD4 24060590 */  addiu   $a2, $zero, 0x0590         ## $a2 = 00000590
/* 01D18 80A21AD8 AC4B0000 */  sw      $t3, 0x0000($v0)           ## 00000000
/* 01D1C 80A21ADC 0C0346A2 */  jal     Matrix_NewMtx              
/* 01D20 80A21AE0 00408825 */  or      $s1, $v0, $zero            ## $s1 = 00000000
/* 01D24 80A21AE4 AE220004 */  sw      $v0, 0x0004($s1)           ## 00000004
/* 01D28 80A21AE8 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 01D2C 80A21AEC 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 01D30 80A21AF0 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 01D34 80A21AF4 AE0C02D0 */  sw      $t4, 0x02D0($s0)           ## 000002D0
/* 01D38 80A21AF8 AC5E0000 */  sw      $s8, 0x0000($v0)           ## 00000000
/* 01D3C 80A21AFC 8ECD0000 */  lw      $t5, 0x0000($s6)           ## 00000000
/* 01D40 80A21B00 01B77021 */  addu    $t6, $t5, $s7              
/* 01D44 80A21B04 01C17821 */  addu    $t7, $t6, $at              
/* 01D48 80A21B08 AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000004
.L80A21B0C:
/* 01D4C 80A21B0C 26730001 */  addiu   $s3, $s3, 0x0001           ## $s3 = 00000001
/* 01D50 80A21B10 00139C00 */  sll     $s3, $s3, 16               
/* 01D54 80A21B14 00139C03 */  sra     $s3, $s3, 16               
/* 01D58 80A21B18 2A610028 */  slti    $at, $s3, 0x0028           
/* 01D5C 80A21B1C 1420FF7F */  bne     $at, $zero, .L80A2191C     
/* 01D60 80A21B20 2652003C */  addiu   $s2, $s2, 0x003C           ## $s2 = 000002B0
/* 01D64 80A21B24 3C0680A2 */  lui     $a2, %hi(D_80A21CB0)       ## $a2 = 80A20000
/* 01D68 80A21B28 24C61CB0 */  addiu   $a2, $a2, %lo(D_80A21CB0)  ## $a2 = 80A21CB0
/* 01D6C 80A21B2C 27A4009C */  addiu   $a0, $sp, 0x009C           ## $a0 = FFFFFFDC
/* 01D70 80A21B30 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 01D74 80A21B34 0C031AD5 */  jal     Graph_CloseDisps              
/* 01D78 80A21B38 24070596 */  addiu   $a3, $zero, 0x0596         ## $a3 = 00000596
/* 01D7C 80A21B3C 8FBF0064 */  lw      $ra, 0x0064($sp)           
/* 01D80 80A21B40 D7B40038 */  ldc1    $f20, 0x0038($sp)          
/* 01D84 80A21B44 8FB00040 */  lw      $s0, 0x0040($sp)           
/* 01D88 80A21B48 8FB10044 */  lw      $s1, 0x0044($sp)           
/* 01D8C 80A21B4C 8FB20048 */  lw      $s2, 0x0048($sp)           
/* 01D90 80A21B50 8FB3004C */  lw      $s3, 0x004C($sp)           
/* 01D94 80A21B54 8FB40050 */  lw      $s4, 0x0050($sp)           
/* 01D98 80A21B58 8FB50054 */  lw      $s5, 0x0054($sp)           
/* 01D9C 80A21B5C 8FB60058 */  lw      $s6, 0x0058($sp)           
/* 01DA0 80A21B60 8FB7005C */  lw      $s7, 0x005C($sp)           
/* 01DA4 80A21B64 8FBE0060 */  lw      $s8, 0x0060($sp)           
/* 01DA8 80A21B68 03E00008 */  jr      $ra                        
/* 01DAC 80A21B6C 27BD00C0 */  addiu   $sp, $sp, 0x00C0           ## $sp = 00000000
