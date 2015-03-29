//------------------------------------------------------------------------------
//
// Элементы подписи с функцией сравнения
// Максим Рябочкин, 28.03.2015
//
//------------------------------------------------------------------------------

(function (root, factory) {

    //
    // Module imports and exports
    // 
    if (typeof define === 'function' && define.amd) {
        define(factory);
    } else if (typeof exports === 'object') {
        module.exports = factory();
    } else {
        root.Signature = factory();
    }

}(this, function () {

    // Построение сплайна
    // x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
    // y - значения функции в узлах сетки
    function buildSpline(x, y) {
        // Количество узлов сетки
        var n = Math.min(x.length, y.length);

        // Инициализация массива сплайнов
        var sp1ines = [];
        for (var i = 0; i < n; ++i)
            sp1ines[i] = {x: x[i], a: y[i], b: 0, c: 0, d: 0};

        // Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
        // Вычисление прогоночных коэффициентов - прямой ход метода прогонки
        var alpha = [], beta = [];
        alpha[0] = beta[0] = 0;
        for (var i = 1; i < n - 1; ++i) {
            var hi = x[i] - x[i - 1],
                    hi1 = x[i + 1] - x[i],
                    A = hi,
                    C = 2 * (hi + hi1),
                    B = hi1,
                    F = 6 * ((y[i + 1] - y[i]) / hi1 - (y[i] - y[i - 1]) / hi),
                    z = (A * alpha[i - 1] + C);
            alpha[i] = -B / z;
            beta[i] = (F - A * beta[i - 1]) / z;
        }

        // Нахождение решения - обратный ход метода прогонки
        for (var i = n - 2; i > 0; --i) {
            sp1ines[i].c = alpha[i] * sp1ines[i + 1].c + beta[i];
        }

        // По известным коэффициентам c[i] находим значения b[i] и d[i]
        for (var i = n - 1; i > 0; --i) {
            var hi = x[i] - x[i - 1];
            sp1ines[i].d = (sp1ines[i].c - sp1ines[i - 1].c) / hi;
            sp1ines[i].b = hi * (2 * sp1ines[i].c + sp1ines[i - 1].c) / 6 + (y[i] - y[i - 1]) / hi;
        }
        return sp1ines;
    }

    // Вычисление значения интерполированной функции в произвольной точке
    function interpolate(sp1ines, x) {

        var n = sp1ines.length, s = {};

        // Если x меньше точки сетки x[0] - пользуемся первым эл-тов массива
        if (x <= sp1ines[0].x) {
            s = sp1ines[0];
            // Если x больше точки сетки x[n - 1] - пользуемся последним эл-том массива
        } else if (x >= sp1ines[n - 1].x) {
            s = sp1ines[n - 1];
            // Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл-та массива
        } else {
            var i = 0;
            var j = n - 1;
            while (i + 1 < j) {
                var k = i + Math.floor((j - i) / 2);
                if (x <= sp1ines[k].x) {
                    j = k;
                } else {
                    i = k;
                }
            }
            s = sp1ines[j];
        }

        var dx = x - s.x;
        // Вычисляем значение сплайна в заданной точке по схеме Горнера 
        return s.a + (s.b + (s.c / 2 + s.d * dx / 6) * dx) * dx;
    }

    // Вычисление коэффициента корреляции 
    // a и b - сравниваемые последовательности
    function correlation(a, b) {
        // Математическое ожидание a[i]
        var a0 = 0, na = a.length;
        for (var i = 0; i < na; i++)
            a0 = a0 + a[i];
        a0 = a0 / na;

        // Математическое ожидание b[i]
        var b0 = 0, nb = b.length;
        for (var i = 0; i < nb; i++)
            b0 = b0 + b[i];
        b0 = b0 / nb;

        // Ковариация a[i] и b[i]
        var cov = 0, da = 0, db = 0;
        for (var i = 0; i < Math.min(na, nb); i++)
            cov = cov + (a[i] - a0) * (b[i] - b0);

        // Квадратичное отклонение (дисперсия) a[i]
        for (var i = 0; i < na; i++)
            da = da + (a[i] - a0) * (a[i] - a0);

        // Квадратичное отклонение (дисперсия) b[i]
        for (var i = 0; i < nb; i++)
            db = db + (b[i] - b0) * (b[i] - b0);

        // Коэффициент корреляции a[i] и b[i]
        return da === 0 && db === 0 ? 1 :
                (da === 0 || db === 0 ? 0 :
                        cov / Math.sqrt(da * db));
    }

    // Преобразование массима объектов в объет массивов
    function toObjectArray(p) {
        var r = {};
        for (var i = 0; i < p.length; i++) {
            var x = p[i];
            for (var name in x) {
                if (!r[name])
                    r[name] = [];
                r[name][i] = x[name];
            }
        }
        return r;
    }

    // Получение массима объектов из объета массивов
    function fromObjectArray(r) {
        var p = [];
        for (var name in r) {
            var h = r[name];
            for (var i = 0; i < h.length; i++) {
                if (typeof h[i] !== 'undefined') {
                    if (!p[i])
                        p[i] = {};
                    p[i][name] = h[i];
                }
            }
        }
        return p;
    }

    // Текущее время
    function now() {
        return (new Date()).getTime();
    }

    // Построение равновременного последовательного ряда координат
    // line - линия
    // size - количество точек в равновременном ряду
    function fillSequence(line, size, sequence) {
        // Преобразование точки в линию по времени
        if (line.length === 1) {
            var point = line[0];
            line = [{x: point.x, y: point.y, t: 0}, {x: point.x, y: point.y, t: 1}];
        }
        // Преобразование в объект массивов
        var a = toObjectArray(line);

        // Приведение временной координаты к интервалу [0..1]
        var n = a.t.length, tmax = a.t[n - 1];
        for (var i = 0; i < n; i++)
            a.t[i] = a.t[i] / tmax;

        // Вычисление сплайнов по каждой координате
        var splices = {
            x: buildSpline(a.t, a.x),
            y: buildSpline(a.t, a.y)
        };

        // Вычисление равновременных последовательностей
        // с границами на концах линий
        for (var i = 0; i <= size; i++) {
            var t = i / size;
            sequence.x.push(interpolate(splices.x, t));
            sequence.y.push(interpolate(splices.y, t));
        }
    }

    // Преобразование в компрессированны формат
    // 
    // Пример кодирования линии
    // [{x:235,y:15,t:0},{x:236,y:13,t:3},{x:241,y:10,t:6},{x:229,y:11,t:11},{x:224,y:10,t:12}]
    // 
    // Промежуточное преобразование к интервальным рядам числ 
    //  x0, dx1, dx2, dx3, dx4; y0, dy1, dy2, dy3, dy4; t0, dt1, dt2, dt3, dt4
    // 235,  +1,  +5, -12,  -5; 15,  -2,  -3,  +1,  -1;  0,  +3,  +3,  +5,  +1
    // 
    // Далее интервальные ряды преобразуются в бинарную кодировку
    // Числа кодируются октетами по 4 бита. Младшие 3 бита определяют число 0..7. 
    // Старший бит - индикатор последнего октета в числе. 1 - последний октет в числе, 0 - читать следующий октет.
    // Значение 0000 определяет смену знака. Начальный знак для каждого ряда '+';
    // Значение 00000000 - окончание ряда.
    // 
    //           235,  +1,  +5,           -12,   -5         ;
    // 0011 0101 1101 1001 1101 0000 0001 1100 1010 0000 0000 
    //       15,       -2,  -3,       +1,        -1         ;  
    // 0001 1111 0000 1010 1011 0000 1001 0000 1001 0000 0000
    //   0,  +3,  +3,  +5,   +1 
    // 1000 1011 1011 1101 1001 
    //
    // Полученная строка преобразуется в base64 или бинарный (binary = true) формат
    //
    function compress(stroke, binary) {
        var bytes = [], byte = 0, m = 0, base = binary ? 8 : 6,
                push2 = function (ch) {
                    byte = (byte << 2) | (ch & 3);
                    m += 2;
                    if (m >= base) {
                        bytes.push(byte);
                        byte = 0;
                        m = 0;
                    }
                },
                push4 = function (ch) {
                    push2(ch >>> 2 & 3);
                    push2(ch & 3);
                },
                pushInt = function (value) {
                    var a = [value & 7];
                    value >>>= 3;
                    while (value > 0) {
                        a.push(value & 7);
                        value >>>= 3;
                    }
                    for (var i = a.length - 1; i >= 0; --i) {
                        if (i === 0)
                            push4(8 | a[i]);
                        else
                            push4(a[i]);
                    }
                },
                pushEnd = function () {
                    if (m > 0) {
                        byte <<= (base - m);
                        bytes.push(byte);
                        byte = 0;
                        m = 0;
                    }
                };

        for (var i = 0, n = stroke.length; i < n; i++) {
            var a = toObjectArray(stroke[i]), h = [a.x, a.y, a.t];
            for (var k = 0, kn = h.length; k < kn; k++) {
                var s = h[k], sign = 1;
                pushInt(s[0]);
                for (var j = 1, jn = s.length; j < jn; j++) {
                    var ds = s[j] - s[j - 1];
                    if (ds !== 0)
                        ds = sign * ds;
                    if (ds < 0) {
                        ds = -ds;
                        sign = -sign;
                        push4(0); // Смена знака
                    }
                    pushInt(ds);
                }
                if (k < kn - 1) {
                    push4(0); // Разделитель рядов
                    push4(0);
                }
            }
            if (i < n - 1) {
                push4(0); // Разделитель рядов
                push4(0);
            }
        }
        pushEnd(); // Окончание последовательности

        // Финальное преобразование
        if (binary) {
            return bytes;
        } else {
            var s = '';
            for (var i = 0; i < bytes.length; i++) {
                var c = bytes[i];
                c = c < 26 ? c + 65 : c < 52 ? c + 71 : c < 62 ? c - 4 :
                        c === 62 ? 43 : c === 63 ? 47 : 65;
                s += String.fromCharCode(c);
            }
            return s;
        }
    }

    // Распаковать данные в компрессированном формате
    function decompress(s) {
        var bytes, base;
        if (typeof s === 'string' || s instanceof String) {
            base = 6;
            bytes = [];
            for (var i = 0; i < s.length; i++) {
                var c = s.charCodeAt(i);
                c = c > 64 && c < 91 ?
                        c - 65 : c > 96 && c < 123 ?
                        c - 71 : c > 47 && c < 58 ?
                        c + 4 : c === 43 ?
                        62 : c === 47 ?
                        63 : 0;
                bytes.push(c);
            }
        } else {
            base = 8;
            bytes = s;
        }
        var byte, m = 0, index = 0, sign = 1,
                read2 = function () {
                    if (m === 0) {
                        if (index >= bytes.length)
                            return 0; // Дополнение нулями
                        byte = bytes[index];
                        index++;
                        m = base;
                    }
                    m -= 2;
                    return byte >>> m & 3;
                },
                read4 = function () {
                    var r = read2() << 2;
                    r |= read2();
                    return r;
                },
                readInt = function () {
                    var value = 0, ch = read4(), i = 0;
                    if (ch === 0) {
                        ch = read4();
                        if (ch === 0)
                            return null; // завершение ряда
                        else
                            sign = -sign;
                    }
                    while ((ch & 8) === 0 && i < 10) { // до последнего октета
                        value = (value << 3) | ch;
                        ch = read4();
                        i++; // выход по неверным данным, максимальная длина 10 октет
                    }
                    value = (value << 3) | (ch & 7);
                    return value === 0 ? value : sign * value;
                };
        var stroke = [];
        while (index < bytes.length) { // До окончания файла
            var line = {}, k = 0, names = ['x', 'y', 't'];
            while (k < 3) {
                var a = [], h = 0, g = 0;
                sign = 1;
                while (true) { // до завершения ряда
                    var g = readInt();
                    if (g === null) // завершение ряда
                        break;
                    h += g;
                    a.push(h);

                }
                line[names[k]] = a;
                k++;
            }
            stroke.push(fromObjectArray(line));
        }
        return stroke;
    }

    // Сравнение подписей с применением статистических методов
    // На выходе коэффициенты корреляции по координатам X и Y
    // Если подписи идентичны, оба коэффициента должны превышать заданное пороговое значение
    function compare(stroke1, stroke2) {
        var a = {x: [], y: []}, b = {x: [], y: []}, empty = [{x: 0, y: 0, t: 0}];

        // Приведение сравниваемых кривых к равно временному интервалу методом 
        // интерполяции кубическим сплайном отдельно для каждой линии
        // Результаты записываются в отдельные последовательности по координатам x и y
        for (var i = 0; i < Math.max(stroke1.length, stroke2.length); i++) {
            var aline = stroke1[i] || empty, bline = stroke2[i] || empty;

            // Коэффициент интерполяции выбирается так, чтобы в результате 
            // получилось более 10 значений на 1 узел для каждой линии
            var size = Math.max(aline.length, bline.length) * 10;

            // Заполнение равновременных ряда для линии
            fillSequence(aline, size, a);
            fillSequence(bline, size, b);
        }

        // Вычисление коэффициентов корреляций по каждой координате
        return {
            x: correlation(a.x, b.x),
            y: correlation(a.y, b.y)
        };
    }

    // Определить, является ли последняя проведенная линия перечеркивающей
    function checkCrossOut(stroke) {
        // Должно быть по крайней мере две линии
        if (stroke.length < 2)
            return false;

        var line = stroke[stroke.length - 1], n = line.length;
        // Координаты X и Y линии не меняют знак
        if (n > 2) {
            var cx = line[0].x - line[1].x, cy = line[0].y - line[1].y;
            for (var i = 2; i < n; i++) {
                var dx = line[i - 1].x - line[i].x, dy = line[i - 1].y - line[i].y;
                if ((cx > 0 && dx < 0) || (cx < 0 && dx > 0) ||
                        (cy > 0 && dy < 0) || (cy < 0 && dy > 0))
                    return false;
            }
        }

        // Граничные координты X линии максимальны и минимальны на всем изображении
        var max = Math.max(line[0].x, line[n - 1].x),
                min = Math.min(line[0].x, line[n - 1].x);
        for (var i = 0; i < stroke.length - 1; i++) {
            var a = stroke[i];
            for (var j = 0; j < a.length; j++) {
                if (max < a[j].x || min > a[j].x)
                    return false;
            }
        }

        return true;
    }

    function extend(obj) {
        for (var i = 1, n = arguments.length; i < n; i++) {
            var ext = arguments[i];
            for (var name in ext)
                obj[name] = ext[name];
        }
    }

    function getStyle(element, styleProp) {
        if (element.currentStyle)
            return element.currentStyle[styleProp];
        else if (window.getComputedStyle)
            return  document.defaultView.getComputedStyle(element, null).getPropertyValue(styleProp);
    }

    function offsetPosition(element) {
        var offsetLeft = 0, offsetTop = 0;
        do {
            offsetLeft += element.offsetLeft;
            offsetTop += element.offsetTop;
        } while (element = element.offsetParent);
        return {offsetLeft: offsetLeft, offsetTop: offsetTop};
    }


    function Signature(element, options) {

        var stroke = [], startTime, curLine, lastPoint, ctx, canvas, 
                changeTimer, self = this;
        options = options || {};

        function resize() {
            canvas.width = element.clientWidth;
            canvas.height = element.clientHeight;
            // Инициализировать контекст
            ctx = canvas.getContext('2d');
            ctx.strokeStyle = options.color || getStyle(element, 'color') || '#000000'; // Colour of the signature 
            ctx.lineWidth = options.thickness || 2; // Thickness of the stroke
            ctx.lineCap = 'Math.round';
            ctx.lineJoin = 'Math.round';
            redraw();
        }

        function change() {
            if ('value' in element)
                element.value = compress(stroke);
            options.change && options.change.apply(self);
        }

        function stop(event) {
            // Отключить обработчики событий
            document.removeEventListener('mouseup', stop, false);
            document.removeEventListener('touchend', stop, false);
            document.removeEventListener('touchcancel', stop, false);
            document.removeEventListener('mousemove', move, false);
            document.removeEventListener('touchmove', move, false);
            // Нарисовать точку для одного узла
            if (curLine.length === 1) {
                ctx.beginPath();
                ctx.moveTo(lastPoint.x, lastPoint.y);
                ctx.lineTo(lastPoint.x + 1, lastPoint.y + 1);
                ctx.stroke();
            }
            // Стереть подпись, если зачеркнули
            if (checkCrossOut(stroke))
                clear();
            else // Таймуат изменения
                changeTimer = setTimeout(change, 1000);
            // Остановить обработку этого события
            if (event && event.target.id === canvas) {
                event.preventDefault();
                event.stopImmediatePropagation();
            }
        }

        function move(event) {
            if (event && event.target === canvas) {
                // Определить координаты точки
                var point, pos = offsetPosition(event.target);
                if (event.touches) {
                    point = {
                        x: event.touches[0].pageX - pos.offsetLeft,
                        y: event.touches[0].pageY - pos.offsetTop,
                        t: now() - startTime
                    };
                } else {
                    point = {
                        x: event.offsetX,
                        y: event.offsetY,
                        t: now() - startTime
                    };
                }
                // Минимальный интервал должен быть 1мс
                if (point.t === lastPoint.t)
                    point.t = lastPoint.t + 1;
                // Нарисовать линию
                ctx.beginPath();
                ctx.moveTo(lastPoint.x, lastPoint.y);
                ctx.lineTo(point.x, point.y);
                ctx.stroke();
                // Добавить точку в текущую линию
                curLine.push(point);
                lastPoint = point;
                // Остановить обработку этого события
                event.preventDefault();
                event.stopImmediatePropagation();
            } else
                stop(event);
        }

        function start(event) {
            if (event && event.target === canvas) {
                // Определить координаты первой точки
                startTime = now();
                var pos = offsetPosition(event.target);
                if (event.touches) {
                    lastPoint = {
                        x: event.touches[0].pageX - pos.offsetLeft,
                        y: event.touches[0].pageY - pos.offsetTop,
                        t: 0
                    };
                } else {
                    lastPoint = {
                        x: event.offsetX,
                        y: event.offsetY,
                        t: 0
                    };
                }
                // Добавить линию в подпись
                curLine = [lastPoint];
                stroke.push(curLine);
                // Остановить обработку этого события
                event.preventDefault();
                event.stopImmediatePropagation();
                clearTimeout(changeTimer);
                // Включить обработчики событий
                document.addEventListener('mouseup', stop, false);
                document.addEventListener('touchend', stop, false);
                document.addEventListener('touchcancel', stop, false);
                document.addEventListener('mousemove', move, false);
                document.addEventListener('touchmove', move, false);
            }
        }

        // Очистить подпись
        function clear() {
            if (ctx) {
                ctx.clearRect(0, 0, canvas.width, canvas.height);
                stroke = [];
                change();
            }
        }

        // Перерисовать подпись
        function redraw() {
            if (ctx && stroke) {
                ctx.clearRect(0, 0, canvas.width, canvas.height);
                for (var i = 0; i < stroke.length; i++) {
                    var line = stroke[i];
                    ctx.beginPath();
                    for (var j = 0; j < line.length; j++)
                        ctx[j === 0 ? 'moveTo' : 'lineTo'](line[j].x, line[j].y);
                    if (line.length === 1)
                        ctx.lineTo(line[0].x + 1, line[0].y + 1);
                    ctx.stroke();
                }
            }
        }

        function remove() {
            window.removeEventListener('resize', resize);
            document.removeEventListener('mousedown', start, false);
            document.removeEventListener('touchstart', start, false);
            element.removeChild(canvas);
        }

        // Преобразование в SVG
        function toSVG(stroke) {
            var sa = [];
            for (var i = 0; i < stroke.length; i++) {
                var line = stroke[i], sline = [];
                for (var j = 0; j < line.length; j++)
                    sline.push('' + line[j].x + ',' + line[j].y);
                if (line.length === 1)
                    sline.push('' + (line[0].x + 1) + ',' + (line[0].y + 1));
                sa.push('		<polyline points="' + sline.join(' ') + '"/>\n');
            }
            return '<?xml version="1.0"?>\n<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n' +
                    '<svg xmlns="http://www.w3.org/2000/svg" width="' + canvas.width + '" height="' + canvas.height + '">\n' +
                    '	<g fill="none" stroke="' + options.color + '" stroke-width="' + options.thickness + '">\n' +
                    sa.join('') +
                    '	</g>\n' +
                    '</svg>\n';
        }

        canvas = document.createElement('canvas');
        // canvas.setAttribute('width', element.clientWidth);
        // canvas.setAttribute('height', element.clientHeight);
        element.appendChild(canvas);
        resize();

        window.addEventListener('resize', resize);
        document.addEventListener('mousedown', start, false);
        document.addEventListener('touchstart', start, false);

        extend(self, {
            // Очистить подпись
            clear: clear,
            // Сравнить подпись с эталоном
            compare: function (etalon) {
                var stroke2;
                if (etalon.stroke)
                    stroke2 = etalon.stroke;
                else
                    stroke2 = decompress(etalon);
                return compare(stroke, stroke2);
            },
            // Получить элемент
            getElement: function () {
                return element;
            },
            // Загрузить из Base64
            fromBase64: function (s) {
                stroke = decompress(s);
                redraw();
                change();
            },
            // Загрузить из JSON
            fromJSON: function (s) {
                stroke = JSON.parse(s).stroke;
                redraw();
                change();
            },
            // Подпись отсутствует
            isEmpty: function () {
                return stroke.length === 0;
            },
            // Перерисовать подпись 
            redraw: redraw,
            // Удалить подпись
            remove: remove,
            // Преобразовать в Base64
            toBase64: function () {
                return compress(stroke);
            },
            // Преобразование в ImageData
            toImageData: function () {
                return ctx.getImageData(0, 0, canvas.width, canvas.height);
            },
            // Преобразовать в JSON
            toJSON: function () {
                return '{"stroke":' + JSON.stringify(stroke) + '}';
            },
            // Преобразовать в SVG
            toSVG: function () {
                return toSVG(stroke);
            }
        });
    }

    return Signature;

}));
