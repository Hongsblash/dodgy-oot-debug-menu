glabel func_8098544C
/* 0086C 8098544C 3C028016 */  lui     $v0, 0x8016                ## $v0 = 80160000
/* 00870 80985450 2442E660 */  addiu   $v0, $v0, 0xE660           ## $v0 = 8015E660
/* 00874 80985454 904E1415 */  lbu     $t6, 0x1415($v0)           ## 8015FA75
/* 00878 80985458 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 0087C 8098545C 24010004 */  addiu   $at, $zero, 0x0004         ## $at = 00000004
/* 00880 80985460 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00884 80985464 15C10019 */  bne     $t6, $at, .L809854CC       
/* 00888 80985468 AFA40028 */  sw      $a0, 0x0028($sp)           
/* 0088C 8098546C 8C4F1360 */  lw      $t7, 0x1360($v0)           ## 8015F9C0
/* 00890 80985470 3C088098 */  lui     $t0, %hi(D_8098786C)       ## $t0 = 80980000
/* 00894 80985474 24180001 */  addiu   $t8, $zero, 0x0001         ## $t8 = 00000001
/* 00898 80985478 29E10004 */  slti    $at, $t7, 0x0004           
/* 0089C 8098547C 10200013 */  beq     $at, $zero, .L809854CC     
/* 008A0 80985480 2508786C */  addiu   $t0, $t0, %lo(D_8098786C)  ## $t0 = 8098786C
/* 008A4 80985484 8CA31C44 */  lw      $v1, 0x1C44($a1)           ## 00001C44
/* 008A8 80985488 AC980260 */  sw      $t8, 0x0260($a0)           ## 00000260
/* 008AC 8098548C ACA81D68 */  sw      $t0, 0x1D68($a1)           ## 00001D68
/* 008B0 80985490 24090002 */  addiu   $t1, $zero, 0x0002         ## $t1 = 00000002
/* 008B4 80985494 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 008B8 80985498 A0491414 */  sb      $t1, 0x1414($v0)           ## 8015FA74
/* 008BC 8098549C 2405006A */  addiu   $a1, $zero, 0x006A         ## $a1 = 0000006A
/* 008C0 809854A0 0C021344 */  jal     Item_Give              
/* 008C4 809854A4 AFA3001C */  sw      $v1, 0x001C($sp)           
/* 008C8 809854A8 8FAA0028 */  lw      $t2, 0x0028($sp)           
/* 008CC 809854AC 34018000 */  ori     $at, $zero, 0x8000         ## $at = 00008000
/* 008D0 809854B0 8FA3001C */  lw      $v1, 0x001C($sp)           
/* 008D4 809854B4 85420032 */  lh      $v0, 0x0032($t2)           ## 00000032
/* 008D8 809854B8 00411021 */  addu    $v0, $v0, $at              
/* 008DC 809854BC 00021400 */  sll     $v0, $v0, 16               
/* 008E0 809854C0 00021403 */  sra     $v0, $v0, 16               
/* 008E4 809854C4 A46200B6 */  sh      $v0, 0x00B6($v1)           ## 000000B6
/* 008E8 809854C8 A4620032 */  sh      $v0, 0x0032($v1)           ## 00000032
.L809854CC:
/* 008EC 809854CC 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 008F0 809854D0 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 008F4 809854D4 03E00008 */  jr      $ra                        
/* 008F8 809854D8 00000000 */  nop


