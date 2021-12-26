## [1-7 задание](./1_7_test.ipynb)

## 8. Тригонометрия разнообразная [(решение)](./trygonometry.py)

Ограничение времени 	1 секунда

Ограничение памяти 	256Mb

Ввод 	стандартный ввод или input.txt

Вывод 	стандартный вывод или output.txt

Дан массив размера _n_ целых чисел _ai_. Необходимо посчитать количество различных значений выражения _sin(cos(ai))_.

### Формат ввода

В первой строке содержится одно число _n_ (2≤ _n_ ≤10^5) — длина массива.

Во второй строке содержится n целых чисел ai — элементы массива (−10^9≤ _ai_ ≤10^9).

### Формат вывода

Выведите одно целое число - ответ на задачу.

### Пример

**Ввод**
```
3
-1 2 2
```

**Вывод**
```
2
```


## 9. Двойная рекуррента [(неполное решение)](./9_recur.py)

Ограничение времени 	1 секунда

Ограничение памяти 	256Mb

Ввод 	стандартный ввод или input.txt

Вывод 	стандартный вывод или output.txt

Пусть последовательности _ak_ и _bk_ заданы рекуррентными уравнениями ("рекурренатами"):

_ak = a{k−1}x0+b{k−1}x1,_

_bk = a{k−1}y0+b{k−1}y1._

Для всех целых _k_>0.

Найдите значение _an_, так как оно может быть достаточно большим посчитайте его по модулю 1000000007. 

### Формат ввода

В первой строке содержится одно число _t_ (1≤ _t_ ≤10^4) — количество наборов входных данных. Далее следует описание наборов входных данных.

В единственной строке содержатся 7 чисел _n, a0, b0, x0, x1, y0, y1_

(1 ≤ _n,a0,b0,x0,x1,y0,y1_ ≤10^9). 

### Формат вывода

Для каждого набора входных данных выведите _an_. 

### Пример

**Ввод**
```
5
1 1 1 1 1 1 1
2 2 1 1 1 3 1
3 1 1 1 1 1 1
4 2 1 2 1 3 1
6 3 4 2 5 8 3
```

**Вывод**
```
2
10
8
185
1415852
```

## 10. Наименьшее общее кратное [(решение)](./10_lcm.py)

Ограничение времени 	50 секунд

Ограничение памяти 	1Gb

Ввод 	стандартный ввод или input.txt

Вывод 	стандартный вывод или output.txt

Даны числа _n,m_. Найдите число различных значений у выражения HOK(i,n)/i+HOK(j,m)/j, для всех натуральных _i,j_, где HOK — наименьшее общее кратное пары чисел. 

### Формат ввода

В первой строке содержится одно число _t_ (1≤ _t_ ≤10) — количество наборов входных данных. Далее следует описание наборов входных данных.

В единственной строке содержатся 2 числа _n, m_ (1≤ _n,m_ ≤10^10). 

### Формат вывода

Для каждого набора входных данных выведите число различных значений выражения из условия. 

### Пример

**Ввод**
```
4
2 4
30 60
1 1
235 87
```

**Вывод**
```
5
45
1
16
```

# 11. Синтетическая классификация [(решение)](./aboba.ipynb)

Вам даны синтетические данные для решения задачи классификации на 5 классов. Каждый объект имеет 15 признаков. Классы обозначены числами от 0 до 4.

Ссылка на данные: https://disk.yandex.ru/d/5r0CNeYbxlMwKg

Данные состоят из трёх файлов. Файл train.csv состоит из 15 колонок, соответствующих признакам, пронумерованным числами от 0 до 14, а также единственной колонки "label", содержащей целевую переменную. Файл test.csv имеет такую же структуру, но в нём пропущена колонка "label".

Файл sample_submission.csv содержит пример посылки. Первая строка файла посылки игнорируется, а в каждой следующей строке должно быть написано одно число: label соответствующего объекта.

Ваша посылка будет оцениваться по метрике accuracy. Формула для оценки: 10⋅min(1,2⋅s)10⋅min(1,2⋅s), где s — значение метрики accuracy.

# 12. Предсказание осадков [(решение)](./aboba.ipynb)

На основе информации про текущий день нужно определить, будут ли осадки завтра.

Данные: https://disk.yandex.ru/d/RGFX7bXCWP3LTg

Предсказываемая колонка — RainTomorrow.

В обучающих данных (train) есть целевая переменная RainTomorrow, в тестовых данных (test) она отсутствует.

Результатом должен быть .csv файл с одной колонкой, i-я строка в котором содержит ответ на i-ю строку файла test_data.csv

Рекомендуется сохранять этот файл с помощью pandas с опциями header=None и index=False.

Результат оценивается по метрике F1-score. Ваш итоговый балл равен min(10,15⋅s)min(10,15⋅s), где ss — f1-мера.

# 13. Бабочки [(решение)](./13_butterflies.ipynb)

Вам необходимо решить задачу классификации фотографий бабочек на 50 классов. Каждая фотография представляет из себя цветное изображение размера 224×224224×224.

В качестве ответа загрузите csv-файл с единственной колонкой. Первая строчка колонки — её название. Далее в строчках для каждого изображения из тестового датасета дан ответ — класс, к которому относится фотография.

Данные можно скачать по ссылке: https://disk.yandex.ru/d/sTESJdMpLNHBTA

Метки классов для обучающих и валидационных данных даны в файлах [train,valid]_labels.csv.

Решение оценивается метрикой accuracy. Балл за задачу равен 20⋅max(0,s−0.5)20⋅max(0,s−0.5), где ss — значение метрики accuracy.

# 14. Выбор игры [(решение)](./aboba.ipynb)

Студент Петя, который также является стажером в области машинного обучения, решил выбрать себе компьютерную игру на некоторой платформе. Но он решил поступить нестандартно и выбрать игру с помощью машинного обучения. Для этого Петя выгрузил базу отзывов на игры. Чтобы выбрать лучшую игру, нужно эти отзывы классифицировать: понять, какие отзывы положительные, а какие — отрицательные. Помогите Пете это сделать.

Вам дана обучающая выборка данных reviews_train.csv. В ней две колонки:

    - id: уникальный идентификатор отзыва;

    - review: текст отзыва на игру;

    - like: целевая переменная: оценка игры пользователем, написавшим отзыв. 1 - положительная, -1 — отрицательная.

Вам нужно обучить алгоритм машинного обучения, который по тексту отзыва (review) предсказывает значение оценки (like) для отзывов из тестовой выборки reviews_test.csv. В качестве ответа загрузите csv файл с ответами, где для каждого id отзыва из тестового датасета дан ответ: 1 или -1. Пример файла с ответом — sample_submission.csv.

Решение оценивается метрикой F1 (http://bazhenov.me/blog/2012/07/21/classification-performance-evaluation.html). Балл вычисляется по формуле 10s10s, где ss — значение метрики F1.

Данные можно скачать по этой ссылке: https://disk.yandex.ru/d/N4LTf1hFKqkwqA

# 15. Переключение с подкреплением

Правильное планирование своего времени достаточно важная и сложная задача, когда Алиса это поняла, то решила вести дневник своей активности. Девушка выделила 3 активности: учеба, развлечения и отдых. В дневнике Алисы 5 колонок: текущее состояние её дел; вид деятельности который она выбрала; полученная польза от выбранной активности; новое состояние её дел, час спустя.

Ваша задача - помочь Алисе с переключением между этими тремя видами деятельности, используя предоставленный дневник. Каждому текущему состоянию дел нужно назначить вероятность выбора той или иной деятельности, чтобы максимизировать итоговую суммарную пользу. Сумма вероятностей по всем деятельностям одного состояния дел должна равняться 100%100%, и не иметь отрицательных слагаемых.

### Формат ввода

Ссылка на данные: https://disk.yandex.ru/d/Dns6gGJIkNHq-A

Вам предоставляется файл history.csv, в котором следующие колонки: состояние дел, выбранная деятельность, польза, следующее состояние дел. Файл submit.csv содержит образец посылки.

### Формат вывода

Вам нужно отправить единственный submit.csv файл, со следующими колонками: состояние дел, процент для выбора каждой деятельности 0−20−2. Процент выбора каждой деятельности должен быть целым неотрицательным числом. Суммарный процент по всем деятельностям одного состояния дел должен равняться 100100.

### Примечания

Также известно, что во время испытания вашей методики будет проводиться 100100 шагов выбора деятельности, начиная из каждого возможного состояния дел. Из каждого состояния будет запущено 100100 таких эпизодов. Итоговые баллы будут определяться конечной суммарной усредненной пользой по всем этим запускам. Чтобы получить максимальный балл ваша стратегия переключения должна действовать не хуже решения жюри.
