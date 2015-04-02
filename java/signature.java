/*
 * Сравнение подписи с эталоном на сервере
 * Максим Рябочкин, 03.04.2015
 */

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 *
 * @author Максим
 */
@WebServlet(urlPatterns = {"/signature"})
public class signature extends HttpServlet {

    private class Point {
        Point(double x0, double y0) {
            x = x0;
            y = y0;
        }
        double x;
        double y;
    }

    private class PointList {
        double[] x;
        double[] y;
    }

    private class Node {
        Node(double x0, double y0, double t0) {
            x = x0;
            y = y0;
            t = t0;
        }
        double x;
        double y;
        double t;
    }

    private class NodeList {
        double[] x;
        double[] y;
        double[] t;
    }

    private class Spline {
        Spline(double x0, double a0, double b0, double c0, double d0) {
            x = x0;
            a = a0;
            b = b0;
            c = c0;
            d = d0;
        }
        double x;
        double a;
        double b;
        double c;
        double d;
    }
    
    // Преобразование массима объектов в объет массивов
    private NodeList toObjectArray(Node[] p) {
        NodeList r = new NodeList();
        r.x = new double[p.length];
        r.y = new double[p.length];
        r.t = new double[p.length];
        for (int i = 0; i < p.length; i++) {
            r.x[i] = p[i].x;
            r.y[i] = p[i].y;
            r.t[i] = p[i].t;
        }
        return r;
    }

    // Получение массима объектов из объета массивов
    private Node[] fromObjectArray(NodeList r) {
        Node[] p = new Node[r.t.length];
        for (int i = 0; i < r.t.length; i++) {
            p[i] = new Node(r.x[i], r.y[i], r.t[i]);
        }
        return p;
    }

    // Распаковать данные в компрессированном формате
    private class Decompressor {

        private int[] bytes;
        private final int base = 6;
        private int bt = 0, m = 0, index = 0, sign = 1;

        private int read2() {
            if (m == 0) {
                if (index >= bytes.length) {
                    return 0; // Дополнение нулями
                }
                bt = bytes[index];
                index++;
                m = base;
            }
            m -= 2;
            return bt >>> m & 3;
        }

        private int read4() {
            int r = read2() << 2;
            r |= read2();
            return r;
        }

        private int readInt() {
            int value = 0, ch = read4(), i = 0;
            if (ch == 0) {
                ch = read4();
                if (ch == 0) {
                    return Integer.MAX_VALUE; // завершение ряда
                } else {
                    sign = -sign;
                }
            }
            while ((ch & 8) == 0 && i < 10) { // до последнего октета
                value = (value << 3) | ch;
                ch = read4();
                i++; // выход по неверным данным, максимальная длина 10 октет
            }
            value = (value << 3) | (ch & 7);
            return value == 0 ? value : sign * value;
        }

        double[] readList() {
            List<Double> al = new ArrayList();
            int h = 0, g;
            sign = 1;
            while (true) { // до завершения ряда
                g = readInt();
                if (g == Integer.MAX_VALUE) // завершение ряда
                {
                    break;
                }
                h += g;
                al.add((double) h);
            }
            double[] a = new double[al.size()];
            for (int i = 0; i < a.length; i++)
                a[i] = al.get(i);
            return a;
        }

        public Node[][] execute(String s) {
            bytes = new int[s.length()];
            for (int i = 0; i < s.length(); i++) {
                int c = s.codePointAt(i);
                c = c > 64 && c < 91
                        ? c - 65 : c > 96 && c < 123
                                ? c - 71 : c > 47 && c < 58
                                        ? c + 4 : c == 43
                                                ? 62 : c == 47
                                                        ? 63 : 0;
                bytes[i] = c;
            }

            
            List<Node[]> strokeList = new ArrayList();
            while (index < bytes.length) { // До окончания файла
                NodeList line = new NodeList();
                line.x = readList();
                line.y = readList();
                line.t = readList();
                strokeList.add(fromObjectArray(line));
            }
            return strokeList.toArray(new Node[0][]);
        }
    }
    
    // Построение сплайна
    // x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
    // y - значения функции в узлах сетки
    private Spline[] buildSpline(double[] x, double[] y) {
        // Количество узлов сетки
        int n = Math.min(x.length, y.length);

        // Инициализация массива сплайнов
        Spline[] splines = new Spline[n];
        for (int i = 0; i < n; ++i)
            splines[i] = new Spline(x[i], y[i], 0, 0, 0);

        // Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
        // Вычисление прогоночных коэффициентов - прямой ход метода прогонки
        double[] alpha = new double[n - 1], beta = new double[n - 1];
        alpha[0] = beta[0] = 0;
        for (int i = 1; i < n - 1; ++i) {
            double hi = x[i] - x[i - 1],
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
        for (int i = n - 2; i > 0; --i) {
            splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
        }

        // По известным коэффициентам c[i] находим значения b[i] и d[i]
        for (int i = n - 1; i > 0; --i) {
            double hi = x[i] - x[i - 1];
            splines[i].d = (splines[i].c - splines[i - 1].c) / hi;
            splines[i].b = hi * (2 * splines[i].c + splines[i - 1].c) / 6 + (y[i] - y[i - 1]) / hi;
        }
        return splines;
    }

    // Вычисление значения интерполированной функции в произвольной точке
    private double interpolate(Spline[] splines, double x) {

        int n = splines.length; 
        Spline s;

        // Если x меньше точки сетки x[0] - пользуемся первым эл-тов массива
        if (x <= splines[0].x) {
            s = splines[0];
            // Если x больше точки сетки x[n - 1] - пользуемся последним эл-том массива
        } else if (x >= splines[n - 1].x) {
            s = splines[n - 1];
            // Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл-та массива
        } else {
            int i = 0;
            int j = n - 1;
            while (i + 1 < j) {
                int k = i + ((j - i) / 2);
                if (x <= splines[k].x) {
                    j = k;
                } else {
                    i = k;
                }
            }
            s = splines[j];
        }

        double dx = x - s.x;
        // Вычисляем значение сплайна в заданной точке по схеме Горнера 
        return s.a + (s.b + (s.c / 2 + s.d * dx / 6) * dx) * dx;
    }

    // Вычисление коэффициента корреляции 
    // a и b - сравниваемые последовательности
    double correlation(double[] a, double[] b) {
        // Математическое ожидание a[i]
        double a0 = 0; 
        int na = a.length;
        for (int i = 0; i < na; i++)
            a0 = a0 + a[i];
        a0 = a0 / na;

        // Математическое ожидание b[i]
        double b0 = 0, nb = b.length;
        for (int i = 0; i < nb; i++)
            b0 = b0 + b[i];
        b0 = b0 / nb;

        // Ковариация a[i] и b[i]
        double cov = 0, da = 0, db = 0;
        for (int i = 0; i < Math.min(na, nb); i++)
            cov = cov + (a[i] - a0) * (b[i] - b0);

        // Квадратичное отклонение (дисперсия) a[i]
        for (int i = 0; i < na; i++)
            da = da + (a[i] - a0) * (a[i] - a0);

        // Квадратичное отклонение (дисперсия) b[i]
        for (int i = 0; i < nb; i++)
            db = db + (b[i] - b0) * (b[i] - b0);

        // Коэффициент корреляции a[i] и b[i]
        return da == 0 && db == 0 ? 1 :
                (da == 0 || db == 0 ? 0 :
                        cov / Math.sqrt(da * db));
    }
  
    // Построение равновременного последовательного ряда координат
    // line - линия
    // size - количество точек в равновременном ряду
    private void fillSequence(Node[] line, int size, PointList sequence, int ofs) {
        // Преобразование точки в линию по времени
        if (line.length == 1) {
            Node point = line[0];
            line = new Node[]{new Node(point.x, point.y, 0), new Node(point.x, point.y, 1)};
        }
        // Преобразование в объект массивов
        NodeList a = toObjectArray(line);

        // Приведение временной координаты к интервалу [0..1]
        int n = a.t.length; 
        double tmax = a.t[n - 1];
        for (int i = 0; i < n; i++)
            a.t[i] = a.t[i] / tmax;

        // Вычисление сплайнов по каждой координате
        Spline[] splicesX = buildSpline(a.t, a.x);
        Spline[] splicesY = buildSpline(a.t, a.y);

        // Вычисление равновременных последовательностей
        // с границами на концах линий
        for (int i = 0; i <= size; i++) {
            double t = 1.0 * i / size;
            sequence.x[ofs + i] = interpolate(splicesX, t);
            sequence.y[ofs + i] = interpolate(splicesY, t);
        }
    }
    
    // Сравнение подписей с применением статистических методов
    // На выходе коэффициенты корреляции по координатам X и Y
    // Если подписи идентичны, оба коэффициента должны превышать заданное пороговое значение
    private Point compare(Node[][] stroke1, Node[][] stroke2) {


        Node[] empty = new Node[]{new Node(0, 0, 0)};
        int total = 0;
        for (int i = 0; i < Math.max(stroke1.length, stroke2.length); i++) {
            Node[] aline = stroke1.length > i ? stroke1[i] : empty; 
            Node[] bline = stroke2.length > i ? stroke2[i] : empty;
            total += Math.max(aline.length, bline.length) * 10 + 1;
        }
        PointList a = new PointList();
        a.x = new double[total];
        a.y = new double[total];
        PointList b = new PointList();
        b.x = new double[total];
        b.y = new double[total];
        
        // Приведение сравниваемых кривых к равно временному интервалу методом 
        // интерполяции кубическим сплайном отдельно для каждой линии
        // Результаты записываются в отдельные последовательности по координатам x и y
        int ofs = 0;
        for (int i = 0; i < Math.max(stroke1.length, stroke2.length); i++) {
            Node[] aline = stroke1.length > i ? stroke1[i] : empty; 
            Node[] bline = stroke2.length > i ? stroke2[i] : empty;

            // Коэффициент интерполяции выбирается так, чтобы в результате 
            // получилось более 10 значений на 1 узел для каждой линии
            int size = Math.max(aline.length, bline.length) * 10;
            
            // Заполнение равновременных ряда для линии
            fillSequence(aline, size, a, ofs);
            fillSequence(bline, size, b, ofs);
            ofs += size + 1;
        }

        // Вычисление коэффициентов корреляций по каждой координате
        return new Point(correlation(a.x, b.x), correlation(a.y, b.y));
    }
    
    private Node[][] savedSignature;
     
    /**
     * Processes requests for both HTTP <code>GET</code> and <code>POST</code>
     * methods.
     *
     * @param request servlet request
     * @param response servlet response
     * @throws ServletException if a servlet-specific error occurs
     * @throws IOException if an I/O error occurs
     */
    protected void processRequest(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        String result = "Unknown operation";
        String save = request.getParameter("save");
        if (save != null) {
            savedSignature = (new Decompressor()).execute(save); 
            result = "Sample signature saved";
        } 
        String verify = request.getParameter("verify");
        if (verify != null) {
            Node[][] signature = (new Decompressor()).execute(verify);
            if (savedSignature != null) {
                Point cor = compare(signature, savedSignature);
                double c = Math.min(cor.x, cor.y); 
                String comparison = "different";
                if (c > 0.9)
                    comparison = "identical";
                else if (c > 0.6)
                    comparison = "similar";
                result = comparison.substring(0, 1).toUpperCase() + comparison.substring(1);
                result = result + " (" + Math.round(cor.x*100) + "%, " + Math.round(cor.y*100) + "%)";
                result = result + "\",\"x\":" + cor.x + ",\"y\":" + cor.y + ",\"comparison\":\"" + comparison;
            } else
                result = "Sample signature not yet saved";
        } 
        response.setContentType("application/json;charset=UTF-8");
        try (PrintWriter out = response.getWriter()) {
            out.println("{\"result\":\"" + result + "\"}");
        }
    }

    // <editor-fold defaultstate="collapsed" desc="HttpServlet methods. Click on the + sign on the left to edit the code.">
    /**
     * Handles the HTTP <code>GET</code> method.
     *
     * @param request servlet request
     * @param response servlet response
     * @throws ServletException if a servlet-specific error occurs
     * @throws IOException if an I/O error occurs
     */
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        processRequest(request, response);
    }

    /**
     * Handles the HTTP <code>POST</code> method.
     *
     * @param request servlet request
     * @param response servlet response
     * @throws ServletException if a servlet-specific error occurs
     * @throws IOException if an I/O error occurs
     */
    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        processRequest(request, response);
    }

    /**
     * Returns a short description of the servlet.
     *
     * @return a String containing servlet description
     */
    @Override
    public String getServletInfo() {
        return "Short description";
    }// </editor-fold>

}
