.rdata
glabel D_8094ABE8
    .asciz "../z_boss_tw.c"
    .balign 4

glabel D_8094ABF8
    .asciz "../z_boss_tw.c"
    .balign 4

glabel D_8094AC08
    .asciz "../z_boss_tw.c"
    .balign 4

glabel D_8094AC18
    .asciz "../z_boss_tw.c"
    .balign 4

.text
glabel func_80943028
/* 0A358 80943028 27BDFFA8 */  addiu   $sp, $sp, 0xFFA8           ## $sp = FFFFFFA8
/* 0A35C 8094302C AFBF001C */  sw      $ra, 0x001C($sp)           
/* 0A360 80943030 AFB10018 */  sw      $s1, 0x0018($sp)           
/* 0A364 80943034 AFB00014 */  sw      $s0, 0x0014($sp)           
/* 0A368 80943038 AFA5005C */  sw      $a1, 0x005C($sp)           
/* 0A36C 8094303C 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 0A370 80943040 00808825 */  or      $s1, $a0, $zero            ## $s1 = 00000000
/* 0A374 80943044 3C068095 */  lui     $a2, %hi(D_8094ABE8)       ## $a2 = 80950000
/* 0A378 80943048 24C6ABE8 */  addiu   $a2, $a2, %lo(D_8094ABE8)  ## $a2 = 8094ABE8
/* 0A37C 8094304C 27A40040 */  addiu   $a0, $sp, 0x0040           ## $a0 = FFFFFFE8
/* 0A380 80943050 24071AE5 */  addiu   $a3, $zero, 0x1AE5         ## $a3 = 00001AE5
/* 0A384 80943054 0C031AB1 */  jal     Graph_OpenDisps              
/* 0A388 80943058 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 0A38C 8094305C 0C034213 */  jal     Matrix_Push              
/* 0A390 80943060 00000000 */  nop
/* 0A394 80943064 3C014264 */  lui     $at, 0x4264                ## $at = 42640000
/* 0A398 80943068 44813000 */  mtc1    $at, $f6                   ## $f6 = 57.00
/* 0A39C 8094306C C6240028 */  lwc1    $f4, 0x0028($s1)           ## 00000028
/* 0A3A0 80943070 C62C0024 */  lwc1    $f12, 0x0024($s1)          ## 00000024
/* 0A3A4 80943074 8E26002C */  lw      $a2, 0x002C($s1)           ## 0000002C
/* 0A3A8 80943078 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 0A3AC 8094307C 0C034261 */  jal     Matrix_Translate              
/* 0A3B0 80943080 46062380 */  add.s   $f14, $f4, $f6             
/* 0A3B4 80943084 C62C01C8 */  lwc1    $f12, 0x01C8($s1)          ## 000001C8
/* 0A3B8 80943088 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0A3BC 8094308C 44066000 */  mfc1    $a2, $f12                  
/* 0A3C0 80943090 0C0342A3 */  jal     Matrix_Scale              
/* 0A3C4 80943094 46006386 */  mov.s   $f14, $f12                 
/* 0A3C8 80943098 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0A3CC 8094309C 3C18FA00 */  lui     $t8, 0xFA00                ## $t8 = FA000000
/* 0A3D0 809430A0 2419FFFF */  addiu   $t9, $zero, 0xFFFF         ## $t9 = FFFFFFFF
/* 0A3D4 809430A4 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 0A3D8 809430A8 AE0F02D0 */  sw      $t7, 0x02D0($s0)           ## 000002D0
/* 0A3DC 809430AC AC580000 */  sw      $t8, 0x0000($v0)           ## 00000000
/* 0A3E0 809430B0 AC590004 */  sw      $t9, 0x0004($v0)           ## 00000004
/* 0A3E4 809430B4 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0A3E8 809430B8 3C09DA38 */  lui     $t1, 0xDA38                ## $t1 = DA380000
/* 0A3EC 809430BC 35290003 */  ori     $t1, $t1, 0x0003           ## $t1 = DA380003
/* 0A3F0 809430C0 24480008 */  addiu   $t0, $v0, 0x0008           ## $t0 = 00000008
/* 0A3F4 809430C4 AE0802D0 */  sw      $t0, 0x02D0($s0)           ## 000002D0
/* 0A3F8 809430C8 AC490000 */  sw      $t1, 0x0000($v0)           ## 00000000
/* 0A3FC 809430CC 8FAA005C */  lw      $t2, 0x005C($sp)           
/* 0A400 809430D0 3C058095 */  lui     $a1, %hi(D_8094ABF8)       ## $a1 = 80950000
/* 0A404 809430D4 24A5ABF8 */  addiu   $a1, $a1, %lo(D_8094ABF8)  ## $a1 = 8094ABF8
/* 0A408 809430D8 8D440000 */  lw      $a0, 0x0000($t2)           ## 00000000
/* 0A40C 809430DC 24061AFC */  addiu   $a2, $zero, 0x1AFC         ## $a2 = 00001AFC
/* 0A410 809430E0 0C0346A2 */  jal     Matrix_NewMtx              
/* 0A414 809430E4 AFA20038 */  sw      $v0, 0x0038($sp)           
/* 0A418 809430E8 8FA30038 */  lw      $v1, 0x0038($sp)           
/* 0A41C 809430EC 3C040602 */  lui     $a0, 0x0602                ## $a0 = 06020000
/* 0A420 809430F0 2484F608 */  addiu   $a0, $a0, 0xF608           ## $a0 = 0601F608
/* 0A424 809430F4 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 0A428 809430F8 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0A42C 809430FC 00046900 */  sll     $t5, $a0,  4               
/* 0A430 80943100 000D7702 */  srl     $t6, $t5, 28               
/* 0A434 80943104 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 0A438 80943108 AE0B02D0 */  sw      $t3, 0x02D0($s0)           ## 000002D0
/* 0A43C 8094310C 000E7880 */  sll     $t7, $t6,  2               
/* 0A440 80943110 3C0CDE00 */  lui     $t4, 0xDE00                ## $t4 = DE000000
/* 0A444 80943114 3C188016 */  lui     $t8, 0x8016                ## $t8 = 80160000
/* 0A448 80943118 030FC021 */  addu    $t8, $t8, $t7              
/* 0A44C 8094311C 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 0A450 80943120 AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 0A454 80943124 8F186FA8 */  lw      $t8, 0x6FA8($t8)           ## 80166FA8
/* 0A458 80943128 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 0A45C 8094312C 0081C824 */  and     $t9, $a0, $at              
/* 0A460 80943130 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 0A464 80943134 03194021 */  addu    $t0, $t8, $t9              
/* 0A468 80943138 01014821 */  addu    $t1, $t0, $at              
/* 0A46C 8094313C AC490004 */  sw      $t1, 0x0004($v0)           ## 00000004
/* 0A470 80943140 8FAA005C */  lw      $t2, 0x005C($sp)           
/* 0A474 80943144 0C025011 */  jal     func_80094044              
/* 0A478 80943148 8D440000 */  lw      $a0, 0x0000($t2)           ## 00000000
/* 0A47C 8094314C 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0A480 80943150 3C0CFA00 */  lui     $t4, 0xFA00                ## $t4 = FA000000
/* 0A484 80943154 240D00C8 */  addiu   $t5, $zero, 0x00C8         ## $t5 = 000000C8
/* 0A488 80943158 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 0A48C 8094315C AE0B02D0 */  sw      $t3, 0x02D0($s0)           ## 000002D0
/* 0A490 80943160 3C014370 */  lui     $at, 0x4370                ## $at = 43700000
/* 0A494 80943164 44817000 */  mtc1    $at, $f14                  ## $f14 = 240.00
/* 0A498 80943168 AC4D0004 */  sw      $t5, 0x0004($v0)           ## 00000004
/* 0A49C 8094316C AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 0A4A0 80943170 8E26002C */  lw      $a2, 0x002C($s1)           ## 0000002C
/* 0A4A4 80943174 C62C0024 */  lwc1    $f12, 0x0024($s1)          ## 00000024
/* 0A4A8 80943178 0C034261 */  jal     Matrix_Translate              
/* 0A4AC 8094317C 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 0A4B0 80943180 3C01457A */  lui     $at, 0x457A                ## $at = 457A0000
/* 0A4B4 80943184 44815000 */  mtc1    $at, $f10                  ## $f10 = 4000.00
/* 0A4B8 80943188 C6280050 */  lwc1    $f8, 0x0050($s1)           ## 00000050
/* 0A4BC 8094318C 3C0142C8 */  lui     $at, 0x42C8                ## $at = 42C80000
/* 0A4C0 80943190 44819000 */  mtc1    $at, $f18                  ## $f18 = 100.00
/* 0A4C4 80943194 460A4402 */  mul.s   $f16, $f8, $f10            
/* 0A4C8 80943198 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 0A4CC 8094319C 44817000 */  mtc1    $at, $f14                  ## $f14 = 1.00
/* 0A4D0 809431A0 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0A4D4 809431A4 46128303 */  div.s   $f12, $f16, $f18           
/* 0A4D8 809431A8 44066000 */  mfc1    $a2, $f12                  
/* 0A4DC 809431AC 0C0342A3 */  jal     Matrix_Scale              
/* 0A4E0 809431B0 00000000 */  nop
/* 0A4E4 809431B4 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0A4E8 809431B8 3C0FDA38 */  lui     $t7, 0xDA38                ## $t7 = DA380000
/* 0A4EC 809431BC 35EF0003 */  ori     $t7, $t7, 0x0003           ## $t7 = DA380003
/* 0A4F0 809431C0 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 0A4F4 809431C4 AE0E02D0 */  sw      $t6, 0x02D0($s0)           ## 000002D0
/* 0A4F8 809431C8 AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
/* 0A4FC 809431CC 8FB8005C */  lw      $t8, 0x005C($sp)           
/* 0A500 809431D0 3C058095 */  lui     $a1, %hi(D_8094AC08)       ## $a1 = 80950000
/* 0A504 809431D4 24A5AC08 */  addiu   $a1, $a1, %lo(D_8094AC08)  ## $a1 = 8094AC08
/* 0A508 809431D8 8F040000 */  lw      $a0, 0x0000($t8)           ## 00000000
/* 0A50C 809431DC 24061B0E */  addiu   $a2, $zero, 0x1B0E         ## $a2 = 00001B0E
/* 0A510 809431E0 0C0346A2 */  jal     Matrix_NewMtx              
/* 0A514 809431E4 AFA2002C */  sw      $v0, 0x002C($sp)           
/* 0A518 809431E8 8FA3002C */  lw      $v1, 0x002C($sp)           
/* 0A51C 809431EC 3C040405 */  lui     $a0, 0x0405                ## $a0 = 04050000
/* 0A520 809431F0 24849210 */  addiu   $a0, $a0, 0x9210           ## $a0 = 04049210
/* 0A524 809431F4 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 0A528 809431F8 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0A52C 809431FC 00044900 */  sll     $t1, $a0,  4               
/* 0A530 80943200 00095702 */  srl     $t2, $t1, 28               
/* 0A534 80943204 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 0A538 80943208 AE1902D0 */  sw      $t9, 0x02D0($s0)           ## 000002D0
/* 0A53C 8094320C 000A5880 */  sll     $t3, $t2,  2               
/* 0A540 80943210 3C08DE00 */  lui     $t0, 0xDE00                ## $t0 = DE000000
/* 0A544 80943214 3C0C8016 */  lui     $t4, 0x8016                ## $t4 = 80160000
/* 0A548 80943218 018B6021 */  addu    $t4, $t4, $t3              
/* 0A54C 8094321C 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 0A550 80943220 AC480000 */  sw      $t0, 0x0000($v0)           ## 00000000
/* 0A554 80943224 8D8C6FA8 */  lw      $t4, 0x6FA8($t4)           ## 80166FA8
/* 0A558 80943228 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 0A55C 8094322C 00816824 */  and     $t5, $a0, $at              
/* 0A560 80943230 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 0A564 80943234 018D7021 */  addu    $t6, $t4, $t5              
/* 0A568 80943238 01C17821 */  addu    $t7, $t6, $at              
/* 0A56C 8094323C 0C034221 */  jal     Matrix_Pull              
/* 0A570 80943240 AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000004
/* 0A574 80943244 8FB8005C */  lw      $t8, 0x005C($sp)           
/* 0A578 80943248 3C068095 */  lui     $a2, %hi(D_8094AC18)       ## $a2 = 80950000
/* 0A57C 8094324C 24C6AC18 */  addiu   $a2, $a2, %lo(D_8094AC18)  ## $a2 = 8094AC18
/* 0A580 80943250 27A40040 */  addiu   $a0, $sp, 0x0040           ## $a0 = FFFFFFE8
/* 0A584 80943254 24071B15 */  addiu   $a3, $zero, 0x1B15         ## $a3 = 00001B15
/* 0A588 80943258 0C031AD5 */  jal     Graph_CloseDisps              
/* 0A58C 8094325C 8F050000 */  lw      $a1, 0x0000($t8)           ## 00000000
/* 0A590 80943260 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 0A594 80943264 8FB00014 */  lw      $s0, 0x0014($sp)           
/* 0A598 80943268 8FB10018 */  lw      $s1, 0x0018($sp)           
/* 0A59C 8094326C 03E00008 */  jr      $ra                        
/* 0A5A0 80943270 27BD0058 */  addiu   $sp, $sp, 0x0058           ## $sp = 00000000
