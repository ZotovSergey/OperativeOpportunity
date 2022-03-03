# Текстовое обозначение единиц измерения: секунд
SECONDS_STR_NAME = 'seconds'
# Текстовое обозначение единиц измерения: минут
MINUTES_STR_NAME = 'minutes'
# Текстовое обозначение единиц измерения: часов
HOURS_STR_NAME = 'hours'
# Текстовое обозначение единиц измерения: дней
DAYS_STR_NAME = 'days'
# Текстовое обозначение единиц измерения: недель
WEEKS_STR_NAME = 'weeks'
# Текстовое обозначение единиц измерения: месяцев (30 дней)
MONTHS_STR_NAME = 'months'
# Текстовое обозначение единиц измерения: лет
YEARS_STR_NAME = 'years'
# Текстовое обозначение единиц измерения: десятилетий
TEN_YEARS_STR_NAME = 'ten years'

# Текстовое символьное обозначение единиц измерения: секунд
SECONDS_STR_SYMBOL = 'sec'
# Текстовое символьное обозначение единиц измерения: минут
MINUTES_STR_SYMBOL = 'min'
# Текстовое символьное обозначение единиц измерения: часов
HOURS_STR_SYMBOL = 'h'
# Текстовое символьное обозначение единиц измерения: дней
DAYS_STR_SYMBOL = 'd'
# Текстовое символьное обозначение единиц измерения: недель
WEEKS_STR_SYMBOL = 'w'
# Текстовое символьное обозначение единиц измерения: месяцев (30 дней)
MONTHS_STR_SYMBOL = 'mon'
# Текстовое символьное обозначение единиц измерения: лет
YEARS_STR_SYMBOL = 'y'
# Текстовое символьное обозначение единиц измерения: десятилетий
TEN_YEARS_STR_SYMBOL = '10y'

# Количество секунд в минуте
SECONDS_IN_MINUTES = 60
# Количество секунд в часе
SECONDS_IN_HOURS = 60 * SECONDS_IN_MINUTES
# Количество секунд в дне
SECONDS_IN_DAYS = 24 * SECONDS_IN_HOURS
# Количество секунд в неделе
SECONDS_IN_WEEKS = 7 * SECONDS_IN_DAYS
# Количество секунд в месяце (в 30 днях)
SECONDS_IN_MONTHS = 30 * SECONDS_IN_DAYS
# Количество секунд в году
SECONDS_IN_YEARS = 365 * SECONDS_IN_DAYS

# Методы в этом файле написаны для работы с единицами измерения времени (unit), заданными словами:
#  'seconds' - секунды;
#  'minutes' - минуты;
#  'hours' - часы;
#  'weeks' - недели;
#  'days' - дни;
#  'months' - месяцы;
#  'years' - годы;
#  'ten years' - десятилетия.


def unit_in_symbol(unit):
    """
    @Описание:
        Определяет символ для заданной единицы измерения.
    :param unit: единица измерения времени (String)
    :return:
    """
    if unit == SECONDS_STR_NAME:
        symbol = SECONDS_STR_SYMBOL
    elif unit == MINUTES_STR_NAME:
        symbol = MINUTES_STR_SYMBOL
    elif unit == HOURS_STR_NAME:
        symbol = HOURS_STR_SYMBOL
    elif unit == DAYS_STR_NAME:
        symbol = DAYS_STR_SYMBOL
    elif unit == WEEKS_STR_NAME:
        symbol = WEEKS_STR_SYMBOL
    elif unit == MONTHS_STR_NAME:
        symbol = MONTHS_STR_SYMBOL
    elif unit == YEARS_STR_NAME:
        symbol = YEARS_STR_SYMBOL
    elif unit == TEN_YEARS_STR_SYMBOL:
        symbol = TEN_YEARS_STR_NAME
    else:
        symbol = 'error'

    return symbol


def seconds_to_unit(seconds, unit):
    """
    @Описание:
        Метод переводит секунды в заданную единицу измерения времени.
    :param seconds: количество секунд, которые переводятся в единицы измерения времени unit.
    :param unit: единица измерения времени, в которые переводятся секунды seconds (String).
    :return: количество единиц измерения времени unit, в которые переводятся секунды seconds.
    """
    if unit == SECONDS_STR_NAME:
        seconds = seconds
    elif unit == MINUTES_STR_NAME:
        seconds = seconds_to_minutes(seconds)
    elif unit == HOURS_STR_NAME:
        seconds = seconds_to_hours(seconds)
    elif unit == DAYS_STR_NAME:
        seconds = seconds_to_days(seconds)
    elif unit == WEEKS_STR_NAME:
        seconds = seconds_to_weeks(seconds)
    elif unit == MONTHS_STR_NAME:
        seconds = seconds_to_months(seconds)
    elif unit == YEARS_STR_NAME:
        seconds = seconds_to_years(seconds)
    elif unit == TEN_YEARS_STR_SYMBOL:
        seconds = seconds_to_ten_years(seconds)

    return seconds


def to_get_unit_in_seconds(unit):
    """
    @Описание:
        Метод возвращает количество секунд в одной заданной единице измерения времени.
    :param unit: единица измерения времени, для которой возвращается количество секунд (String).
    :return: количество секунд в одной единице измерения времени unit.
    """
    if unit == SECONDS_STR_NAME:
        seconds_quantity = 1
    elif unit == MINUTES_STR_NAME:
        seconds_quantity = SECONDS_IN_MINUTES
    elif unit == HOURS_STR_NAME:
        seconds_quantity = SECONDS_IN_HOURS
    elif unit == DAYS_STR_NAME:
        seconds_quantity = SECONDS_IN_DAYS
    elif unit == WEEKS_STR_NAME:
        seconds_quantity = SECONDS_IN_WEEKS
    elif unit == MONTHS_STR_NAME:
        seconds_quantity = SECONDS_IN_MONTHS
    elif unit == YEARS_STR_NAME:
        seconds_quantity = SECONDS_IN_YEARS
    elif unit == TEN_YEARS_STR_SYMBOL:
        seconds_quantity = 10 * SECONDS_IN_YEARS
    else:
        seconds_quantity = 1

    return seconds_quantity


def seconds_to_minutes(seconds):
    """
    @Описание:
        Метод переводит заданное количество секунд в минуты.
    :param seconds: количество секунд, которое переводится в минуты.
    :return: заданное количество секунд seconds в минутах.
    """
    return seconds / SECONDS_IN_MINUTES


def seconds_to_hours(seconds):
    """
    @Описание:
        Метод переводит заданное количество секунд в часы.
    :param seconds: количество секунд, которое переводится в часы.
    :return: заданное количество секунд seconds в часах.
    """
    return seconds / SECONDS_IN_HOURS


def seconds_to_days(seconds):
    """
    @Описание:
        Метод переводит заданное количество секунд в дни.
    :param seconds: количество секунд, которое переводится в дни.
    :return: заданное количество секунд seconds в днях.
    """
    return seconds / SECONDS_IN_DAYS


def seconds_to_weeks(seconds):
    """
    @Описание:
        Метод переводит заданное количество секунд в недели.
    :param seconds: количество секунд, которое переводится в недели.
    :return: заданное количество секунд seconds в неделях.
    """
    return seconds / SECONDS_IN_WEEKS


def seconds_to_months(seconds):
    """
    @Описание:
        Метод переводит заданное количество секунд в месяцы.
    :param seconds: количество секунд, которое переводится в месяцы.
    :return: заданное количество секунд seconds в месяцах.
    """
    return seconds / SECONDS_IN_MONTHS


def seconds_to_years(seconds):
    """
    @Описание:
        Метод переводит заданное количество секунд в годы.
    :param seconds: количество секунд, которое переводится в годы.
    :return: заданное количество секунд seconds в годах.
    """
    return seconds / SECONDS_IN_YEARS


def seconds_to_ten_years(seconds):
    """
    @Описание:
        Метод переводит заданное количество секунд в десятилетия.
    :param seconds: количество секунд, которое переводится в десятилетия.
    :return: заданное количество секунд seconds в десятилетиях.
    """
    return seconds / (10 * SECONDS_IN_YEARS)
