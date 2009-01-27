      IMPLICIT          NONE
      REAL              REF_SECOND      / 0 /
      INTEGER           UTOPEN
      INTEGER           UTDEC
      INTEGER           UTICALTIME
      INTEGER           STATUS
      INTEGER           REF_YEAR        / 1990 /
      INTEGER           REF_MONTH       / 1 /
      INTEGER           REF_DAY         / 1 /
      INTEGER           REF_HOUR        / 1 /
      INTEGER           REF_MINUTE      / 0 /
      INTEGER           UTMAKE
      INTEGER           TIMECENTERS_UNIT
      DOUBLEPRECISION   REF_VALUE

      STATUS = UTOPEN('udunits.dat')
      IF (STATUS .NE. 0) THEN
          PRINT *, 'Couldn''t open database: status =', STATUS
          CALL ABORT
      ENDIF

      TIMECENTERS_UNIT = UTMAKE()

      STATUS = UTDEC('2 minutes since 1990-1-1', TIMECENTERS_UNIT)
      IF (STATUS .NE. 0)
     *THEN
          PRINT *, 'UTDEC() =', STATUS
      ELSE
C
C         Reference time is start time plus 1 hour.
C
          STATUS = UTICALTIME(REF_YEAR, REF_MONTH, REF_DAY, REF_HOUR, 
     *                    REF_MINUTE, REF_SECOND, TIMECENTERS_UNIT,
     *                    REF_VALUE)
C
C         Number of time intervals between start and reference times:
C
          IF (STATUS .NE. 0) THEN
              PRINT *, 'UTICALTIME() =', STATUS
          ELSE
              IF (30 .NE. REF_VALUE) THEN
                  PRINT *, 'Incorrect result:', REF_VALUE
              ELSE
                  PRINT *, 'Correct result'
              ENDIF
          ENDIF
      ENDIF
      END
