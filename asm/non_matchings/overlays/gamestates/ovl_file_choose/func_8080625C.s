glabel func_8080625C
/* 0251C 8080625C 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 02520 80806260 3C030001 */  lui     $v1, 0x0001                ## $v1 = 00010000
/* 02524 80806264 34638000 */  ori     $v1, $v1, 0x8000           ## $v1 = 00018000
/* 02528 80806268 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 0252C 8080626C AFB00020 */  sw      $s0, 0x0020($sp)           
/* 02530 80806270 00834021 */  addu    $t0, $a0, $v1              
/* 02534 80806274 850E4A38 */  lh      $t6, 0x4A38($t0)           ## 00004A38
/* 02538 80806278 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0253C 8080627C 51C00006 */  beql    $t6, $zero, .L80806298     
/* 02540 80806280 96020020 */  lhu     $v0, 0x0020($s0)           ## 00000020
/* 02544 80806284 948F0020 */  lhu     $t7, 0x0020($a0)           ## 00000020
/* 02548 80806288 31F89000 */  andi    $t8, $t7, 0x9000           ## $t8 = 00000000
/* 0254C 8080628C 57000007 */  bnel    $t8, $zero, .L808062AC     
/* 02550 80806290 85094A46 */  lh      $t1, 0x4A46($t0)           ## 00004A46
/* 02554 80806294 96020020 */  lhu     $v0, 0x0020($s0)           ## 00000020
.L80806298:
/* 02558 80806298 2401BFFF */  addiu   $at, $zero, 0xBFFF         ## $at = FFFFBFFF
/* 0255C 8080629C 0041C827 */  nor     $t9, $v0, $at              
/* 02560 808062A0 1720001E */  bne     $t9, $zero, .L8080631C     
/* 02564 808062A4 304D9000 */  andi    $t5, $v0, 0x9000           ## $t5 = 00000000
/* 02568 808062A8 85094A46 */  lh      $t1, 0x4A46($t0)           ## 00004A46
.L808062AC:
/* 0256C 808062AC 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 02570 808062B0 00300821 */  addu    $at, $at, $s0              
/* 02574 808062B4 A429CA38 */  sh      $t1, -0x35C8($at)          ## 0001CA38
/* 02578 808062B8 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 0257C 808062BC 00300821 */  addu    $at, $at, $s0              
/* 02580 808062C0 240A0006 */  addiu   $t2, $zero, 0x0006         ## $t2 = 00000006
/* 02584 808062C4 A42ACA66 */  sh      $t2, -0x359A($at)          ## 0001CA66
/* 02588 808062C8 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 0258C 808062CC 00300821 */  addu    $at, $at, $s0              
/* 02590 808062D0 240B0019 */  addiu   $t3, $zero, 0x0019         ## $t3 = 00000019
/* 02594 808062D4 A42BCA3E */  sh      $t3, -0x35C2($at)          ## 0001CA3E
/* 02598 808062D8 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 0259C 808062DC 24020008 */  addiu   $v0, $zero, 0x0008         ## $v0 = 00000008
/* 025A0 808062E0 3C078013 */  lui     $a3, 0x8013                ## $a3 = 80130000
/* 025A4 808062E4 00300821 */  addu    $at, $at, $s0              
/* 025A8 808062E8 3C0C8013 */  lui     $t4, 0x8013                ## $t4 = 80130000
/* 025AC 808062EC 24E733E0 */  addiu   $a3, $a3, 0x33E0           ## $a3 = 801333E0
/* 025B0 808062F0 A422CA50 */  sh      $v0, -0x35B0($at)          ## 0001CA50
/* 025B4 808062F4 258C33E8 */  addiu   $t4, $t4, 0x33E8           ## $t4 = 801333E8
/* 025B8 808062F8 3C058013 */  lui     $a1, 0x8013                ## $a1 = 80130000
/* 025BC 808062FC 24A533D4 */  addiu   $a1, $a1, 0x33D4           ## $a1 = 801333D4
/* 025C0 80806300 AFAC0014 */  sw      $t4, 0x0014($sp)           
/* 025C4 80806304 AFA70010 */  sw      $a3, 0x0010($sp)           
/* 025C8 80806308 2404483C */  addiu   $a0, $zero, 0x483C         ## $a0 = 0000483C
/* 025CC 8080630C 0C03DCE3 */  jal     Audio_PlaySoundGeneral
              
/* 025D0 80806310 24060004 */  addiu   $a2, $zero, 0x0004         ## $a2 = 00000004
/* 025D4 80806314 10000047 */  beq     $zero, $zero, .L80806434   
/* 025D8 80806318 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L8080631C:
/* 025DC 8080631C 11A0002B */  beq     $t5, $zero, .L808063CC     
/* 025E0 80806320 3C078013 */  lui     $a3, 0x8013                ## $a3 = 80130000
/* 025E4 80806324 850E4A46 */  lh      $t6, 0x4A46($t0)           ## 00004A46
/* 025E8 80806328 3C0D8013 */  lui     $t5, 0x8013                ## $t5 = 80130000
/* 025EC 8080632C 24E733E0 */  addiu   $a3, $a3, 0x33E0           ## $a3 = 801333E0
/* 025F0 80806330 000E7840 */  sll     $t7, $t6,  1               
/* 025F4 80806334 020FC021 */  addu    $t8, $s0, $t7              
/* 025F8 80806338 0303C821 */  addu    $t9, $t8, $v1              
/* 025FC 8080633C A7204A86 */  sh      $zero, 0x4A86($t9)         ## 00004A86
/* 02600 80806340 85094A46 */  lh      $t1, 0x4A46($t0)           ## 00004A46
/* 02604 80806344 25AD33E8 */  addiu   $t5, $t5, 0x33E8           ## $t5 = 801333E8
/* 02608 80806348 3C058013 */  lui     $a1, 0x8013                ## $a1 = 80130000
/* 0260C 8080634C 00095040 */  sll     $t2, $t1,  1               
/* 02610 80806350 020A5821 */  addu    $t3, $s0, $t2              
/* 02614 80806354 01631021 */  addu    $v0, $t3, $v1              
/* 02618 80806358 844C4A86 */  lh      $t4, 0x4A86($v0)           ## 00004A86
/* 0261C 8080635C 24A533D4 */  addiu   $a1, $a1, 0x33D4           ## $a1 = 801333D4
/* 02620 80806360 240428BA */  addiu   $a0, $zero, 0x28BA         ## $a0 = 000028BA
/* 02624 80806364 A44C4A28 */  sh      $t4, 0x4A28($v0)           ## 00004A28
/* 02628 80806368 AFAD0014 */  sw      $t5, 0x0014($sp)           
/* 0262C 8080636C AFA70010 */  sw      $a3, 0x0010($sp)           
/* 02630 80806370 0C03DCE3 */  jal     Audio_PlaySoundGeneral
              
/* 02634 80806374 24060004 */  addiu   $a2, $zero, 0x0004         ## $a2 = 00000004
/* 02638 80806378 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 0263C 8080637C 24020008 */  addiu   $v0, $zero, 0x0008         ## $v0 = 00000008
/* 02640 80806380 00300821 */  addu    $at, $at, $s0              
/* 02644 80806384 A422CA50 */  sh      $v0, -0x35B0($at)          ## 0001CA50
/* 02648 80806388 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 0264C 8080638C 00300821 */  addu    $at, $at, $s0              
/* 02650 80806390 240E001B */  addiu   $t6, $zero, 0x001B         ## $t6 = 0000001B
/* 02654 80806394 A42ECA3E */  sh      $t6, -0x35C2($at)          ## 0001CA3E
/* 02658 80806398 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 0265C 8080639C 00300821 */  addu    $at, $at, $s0              
/* 02660 808063A0 A422CA66 */  sh      $v0, -0x359A($at)          ## 0001CA66
/* 02664 808063A4 3C014348 */  lui     $at, 0x4348                ## $at = 43480000
/* 02668 808063A8 44816000 */  mtc1    $at, $f12                  ## $f12 = 200.00
/* 0266C 808063AC 240500FF */  addiu   $a1, $zero, 0x00FF         ## $a1 = 000000FF
/* 02670 808063B0 24060014 */  addiu   $a2, $zero, 0x0014         ## $a2 = 00000014
/* 02674 808063B4 0C02A800 */  jal     func_800AA000              
/* 02678 808063B8 24070096 */  addiu   $a3, $zero, 0x0096         ## $a3 = 00000096
/* 0267C 808063BC 240F000F */  addiu   $t7, $zero, 0x000F         ## $t7 = 0000000F
/* 02680 808063C0 3C018081 */  lui     $at, %hi(D_808124A0)       ## $at = 80810000
/* 02684 808063C4 1000001A */  beq     $zero, $zero, .L80806430   
/* 02688 808063C8 A42F24A0 */  sh      $t7, %lo(D_808124A0)($at)  
.L808063CC:
/* 0268C 808063CC 85024ABC */  lh      $v0, 0x4ABC($t0)           ## 00004ABC
/* 02690 808063D0 3C078013 */  lui     $a3, 0x8013                ## $a3 = 80130000
/* 02694 808063D4 24E733E0 */  addiu   $a3, $a3, 0x33E0           ## $a3 = 801333E0
/* 02698 808063D8 04400003 */  bltz    $v0, .L808063E8            
/* 0269C 808063DC 00021823 */  subu    $v1, $zero, $v0            
/* 026A0 808063E0 10000001 */  beq     $zero, $zero, .L808063E8   
/* 026A4 808063E4 00401825 */  or      $v1, $v0, $zero            ## $v1 = 00000000
.L808063E8:
/* 026A8 808063E8 2861001E */  slti    $at, $v1, 0x001E           
/* 026AC 808063EC 14200010 */  bne     $at, $zero, .L80806430     
/* 026B0 808063F0 24044839 */  addiu   $a0, $zero, 0x4839         ## $a0 = 00004839
/* 026B4 808063F4 3C188013 */  lui     $t8, 0x8013                ## $t8 = 80130000
/* 026B8 808063F8 271833E8 */  addiu   $t8, $t8, 0x33E8           ## $t8 = 801333E8
/* 026BC 808063FC 3C058013 */  lui     $a1, 0x8013                ## $a1 = 80130000
/* 026C0 80806400 24A533D4 */  addiu   $a1, $a1, 0x33D4           ## $a1 = 801333D4
/* 026C4 80806404 AFB80014 */  sw      $t8, 0x0014($sp)           
/* 026C8 80806408 24060004 */  addiu   $a2, $zero, 0x0004         ## $a2 = 00000004
/* 026CC 8080640C AFA70010 */  sw      $a3, 0x0010($sp)           
/* 026D0 80806410 0C03DCE3 */  jal     Audio_PlaySoundGeneral
              
/* 026D4 80806414 AFA80028 */  sw      $t0, 0x0028($sp)           
/* 026D8 80806418 8FA80028 */  lw      $t0, 0x0028($sp)           
/* 026DC 8080641C 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 026E0 80806420 00300821 */  addu    $at, $at, $s0              
/* 026E4 80806424 85194A38 */  lh      $t9, 0x4A38($t0)           ## 00004A38
/* 026E8 80806428 3B290001 */  xori    $t1, $t9, 0x0001           ## $t1 = 00000001
/* 026EC 8080642C A429CA38 */  sh      $t1, -0x35C8($at)          ## 0001CA38
.L80806430:
/* 026F0 80806430 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80806434:
/* 026F4 80806434 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 026F8 80806438 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 026FC 8080643C 03E00008 */  jr      $ra                        
/* 02700 80806440 00000000 */  nop


